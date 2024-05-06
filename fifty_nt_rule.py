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
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


from CDS_annotation_present import handle_cds_transcripts
from cds_determination_protein_coding import determine_cds
from apply_50nt import get_length_last_exon, calculate_50nt_rule
from helper_functions_cds_determination import get_fasta_tid
from calculate_sequence_features import get_stop_codon_identity, get_base_after_stop,\
get_GC_content_in30bp_ribo_window, get_number_of_exons_transcript, find_UPF1_motifs_in3prime, \
get_GC_content_in15bp_ribo_window

def main():
    start_time = time.time()
    #read in gtf file
    custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

    #get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    transcript_ids = sorted(transcript_ids)
    print('number of transcripts in custom gtf: ', len(transcript_ids))

    #build dataframe to store the computed features
    NMD_features_df = pd.DataFrame(columns = ['50_nt', 'last_exon_length',
                                                't_length', 'start_ORF', 'end_ORF', 
                                                'exon_with_stop_length',
                                                'stop_TGA', 'stop_TAA'], #stop being TAG is if both others are 0
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
    
    #get transcript and ORF sequences for the transcripts with CDS
    transcript_ids_wo_cds_set = set(transcript_ids_wo_cds)
    tids_with_cds = [tid for tid in transcript_ids if tid not in transcript_ids_wo_cds_set]
    tids_with_cds_set = set(tids_with_cds)
    CDS_seqs = []
    transcript_seqs = []
    if len(tids_with_cds_set) > 0:
        transcripts_with_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                                if k in tids_with_cds_set}
        transcript_seqs = get_fasta_tid(transcripts_with_cds, sys.argv[3], seq_type = 'exon')
        CDS_seqs = get_fasta_tid(transcripts_with_cds, sys.argv[3], seq_type = 'CDS', plus_stop = True)

    start_no_cds = time.time() - start_time
    print('Time for start and handling cds annotated transcripts', start_no_cds)

    #################################################################################
    ###handle transcripts with no CDS annoation: CDS determintation, then 50 nt rule#
    #################################################################################
    if len(transcript_ids_wo_cds) > 0:
        print('transcripts for which no CDS annotation was given in the custom gtf: ',\
               len(transcript_ids_wo_cds))
        output_name = sys.argv[1].split('/')[-1][:-4]
        #determine CDS for source transcripts
        transcripts_calculated_CDS, sequences = determine_cds(transcript_gtftk_object, transcript_ids_wo_cds,\
                sys.argv[2], sys.argv[3], output_name)
        
        #drop all for which no ORF was determined
        transcripts_calculated_CDS = transcripts_calculated_CDS.dropna(subset = 'name')
        
        #get the length of the last exon per transcript for 50 nt rule (from the custom gtf)
        last_exon_length_dict = get_length_last_exon(transcripts_calculated_CDS['tid'].to_list(), custom_gtf)
        transcripts_calculated_CDS['last_exon_length'] = transcripts_calculated_CDS['tid'].map(last_exon_length_dict)
        
        #sequences = SeqIO.parse(f'Output/{output_name}_transcripts_of_interest.fasta', 'fasta')
        #sequences = [seq for seq in sequences]

        #apply 50nt rule
        transcripts_calculated_CDS = calculate_50nt_rule(transcripts_calculated_CDS, sequences)
    
    
        # Display the columns present in each DataFrame
        #print("Columns in NMD_features_df:", NMD_features_df.columns)
        #print("Columns in transcripts_calculated_CDS:", transcripts_calculated_CDS.columns)
        start =  time.time() 
        print(transcripts_calculated_CDS.columns)
        
        print('1', time.time()-start)
        transcripts_calculated_CDS['ORF_id'] = transcripts_calculated_CDS['name'].str.split('|').str[1]
        print('2', time.time()-start)
        transcripts_calculated_CDS['ORF_id'] = transcripts_calculated_CDS['ORF_id'].astype('string')
        print('3', time.time()-start)
        ORFs = SeqIO.parse(f'Output/{output_name}_ORFS_protein_coding_genes.fasta', "fasta")
        ORF_list = []
        ORF_id_set = set(transcripts_calculated_CDS['ORF_id'].unique())
        for ORF in ORFs:
            if ORF.id in ORF_id_set:
                ORF_list.append(ORF)
        ORFs = ORF_list
        #ORFs = [ORF for ORF in ORFs if ORF.id in set(transcripts_calculated_CDS['ORF_id'].unique())]
        #this step seems to be the problem...
        print('4', time.time()-start)
        
        CDS_seqs = CDS_seqs + ORFs


    
        #join with df for the CDS annotated transcripts
        NMD_features_df =pd.concat([transcripts_calculated_CDS, NMD_features_df], axis=0, join='outer')
        NMD_features_df['50_nt'] = np.where(NMD_features_df['50_nt'].isna(),\
                                 NMD_features_df['50_nt_rule'], NMD_features_df['50_nt'])
        NMD_features_df.drop(['50_nt_rule'], axis = 1, inplace = True)
        CDS_merge =  time.time() - start_time - start_no_cds
        print("Time until CDS no CDS merged:", CDS_merge)

        

        transcript_seqs = transcript_seqs + sequences

        seq_time =  time.time() - start_time - start_no_cds - CDS_merge
        print("time needed to get seqs and CDS seqs for feature extraction:", seq_time)
        total_time_CDS =  time.time() - start_time - start_no_cds
        print("total time for CDS with merge dfs and calc 50nt:", total_time_CDS)

    NMD_features_df = NMD_features_df.dropna(subset = ['t_length', 'start_ORF'])
    #calculate features used for the NMD classifier
    NMD_features_df = get_stop_codon_identity(CDS_seqs, NMD_features_df)
    NMD_features_df['3_UTR_length'] = NMD_features_df['t_length']-NMD_features_df['end_ORF']
    NMD_features_df['5_UTR_length'] = NMD_features_df['start_ORF']
    NMD_features_df['distance_stop_from_start'] = NMD_features_df['end_ORF'] - NMD_features_df['start_ORF'] - 2
    NMD_features_df['stop_150bp_from_start'] = np.where(NMD_features_df['distance_stop_from_start'] > 150, 0, 1)
    NMD_features_df = get_base_after_stop(transcript_seqs, NMD_features_df)
    NMD_features_df = get_GC_content_in30bp_ribo_window(transcript_seqs, NMD_features_df)
    NMD_features_df = get_GC_content_in15bp_ribo_window(transcript_seqs, NMD_features_df)
    NMD_features_df = get_number_of_exons_transcript(transcript_seqs, NMD_features_df)
    NMD_features_df = find_UPF1_motifs_in3prime(transcript_seqs, NMD_features_df)
    
    print(NMD_features_df.head())
    print('number of transcripts for which 50 nt rule was calculated: ', sum(NMD_features_df['50_nt'].notna()))
    print('number of transcripts found to be NMD sensitive (no escape) by 50nt rule ', \
          sum(NMD_features_df['50_nt'][NMD_features_df['50_nt'].notna()]))
    
    if len(transcript_ids_wo_cds) > 0:
        NMD_features_df.drop(['gid', 'gid_target', 'ORF_nr', 'target_length', 'name',\
                              'target_coverage_percentage', 'protein_overlap_perc', 'overlap',\
                                'protein_overlap_aa'], axis = 1, inplace = True)
    #write results to csv
    NMD_features_df.to_csv(sys.argv[4])
    return  0





if __name__ == '__main__':
    main()
