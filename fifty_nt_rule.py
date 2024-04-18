#fifty_nt_rule.py determines the CDS for given transcripts in a gtf
#(if not supplied) and applies the 50 nt rule to them
#the output is a csv file where the 50 nt rule is noted for each transcript
#together with the information of the chosen CDS and if applicable
#the corresponding protein cosing target transcript used for the CDS

#uasge: python fifty_nt_rule.py custom.gtf Ensembl.gtf genome.fasta output_name.csv

import sys
import time

import numpy as np
import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF


from CDS_annotation_present import handle_cds_transcripts
from cds_determination_protein_coding import determine_cds
from apply_50nt import get_length_last_exon, calculate_50nt_rule

def main():
    start_time = time.time()
    #read in gtf file
    custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

    #get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    transcript_ids = sorted(transcript_ids)
    print('number of transcripts in custom gtf: ', len(transcript_ids))

    #build dataframe to store the computed features
    NMD_features_df = pd.DataFrame(columns = ['50_nt', 'has_cds', 'last_exon_length',  't_length', 'start_ORF', 'end_ORF'],
                                    index = transcript_ids)

    #get exons, CDS and stop codons from the gtf object as dict
    transcript_gtftk_object = custom_gtf\
        .select_by_key('feature', 'CDS,exon,stop_codon')\
        .extract_data('transcript_id,start,end,exon_number,feature,strand,chrom,gene_id,score',
                       as_dict_of_merged_list=True)
    
    #handling of transcripts with CDS anno in custom gtf
    #not down transcript ids for which no CDS anno found
    transcript_ids_wo_cds, NMD_features_df =\
    handle_cds_transcripts(transcript_gtftk_object, transcript_ids, NMD_features_df)
    print()

    start_no_cds = time.time() - start_time
    print('Time for start and handling cds annotated transcripts', start_no_cds)

    #################################################################################
    ###handle transcripts with no CDS annoation: CDS determintation, then 50 nt rule#
    #################################################################################
    if len(transcript_ids_wo_cds) > 0:
        print('transcripts for which no CDS annotation was given in the custom gtf: ',\
               len(transcript_ids_wo_cds))
        
        #determine CDS for source transcripts
        transcripts_calculated_CDS, sequences = determine_cds(transcript_gtftk_object, transcript_ids_wo_cds,\
                sys.argv[2], sys.argv[3], sys.argv[1].split('/')[-1][:-4])
        
        #get the length of the last exon per transcript for 50 nt rule (from the custom gtf)
        last_exon_length_dict = get_length_last_exon(transcripts_calculated_CDS['tid'].to_list(), custom_gtf)
        transcripts_calculated_CDS['last_exon_length'] = transcripts_calculated_CDS['tid'].map(last_exon_length_dict)
        
        #apply 50nt rule
        transcripts_calculated_CDS = calculate_50nt_rule(transcripts_calculated_CDS, sequences)
    
    
        # Display the columns present in each DataFrame
        print("Columns in NMD_features_df:", NMD_features_df.columns)
        print("Columns in transcripts_calculated_CDS:", transcripts_calculated_CDS.columns)

        #join with df for the CDS annotated transcripts
        NMD_features_df =pd.concat([transcripts_calculated_CDS, NMD_features_df], axis=0, join='outer')
        #NMD_features_df = NMD_features_df.merge(transcripts_calculated_CDS, how='outer', on =['last_exon_length', 'start_ORF',\
         #                                                        'end_ORF', 't_length', 'has_cds', '50_nt'])#, 
        NMD_features_df['50_nt'] = np.where(NMD_features_df['50_nt'].isna(),\
                                 NMD_features_df['50_nt_rule'], NMD_features_df['50_nt'])
        NMD_features_df.drop(['50_nt_rule'], axis = 1, inplace = True)
        
    
    print(NMD_features_df.head())
    print(NMD_features_df.tail())
    print('number of transcripts for which 50 nt rule was calculated: ', sum(NMD_features_df['50_nt'].notna()))
    print('number of transcripts found to be NMD sensitive (no escape) by 50nt rule ', \
          sum(NMD_features_df['50_nt'][NMD_features_df['50_nt'].notna()]))
    
    #write results to csv
    NMD_features_df.to_csv(sys.argv[4])
    return  0





if __name__ == '__main__':
    main()
