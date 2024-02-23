import sys
import time

import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from CDS_annotation_present import handle_cds_transcripts
from cds_determination import determine_cds

def main():
    start_time = time.time()
    #read in gtf files
    #ensembl_gtf = GTF(sys.argv[2], check_ensembl_format=False)
    custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

    #get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    transcript_ids = sorted(transcript_ids)
    print("number of transcripts to investigate", len(transcript_ids))


    #build dataframe to store the computed features
    NMD_features_df = pd.DataFrame(columns = ["50_nt", "has_cds"],
                                    index = transcript_ids)

    #select information on transcripts as a list
    #transcript_string = get_transcript_string(transcript_ids)
    transcript_gtftk_object = custom_gtf\
        .select_by_key('feature', "CDS,exon,stop_codon")\
        .extract_data('transcript_id,start,end,exon_number,feature,strand,chrom,gene_id,score',
                       as_dict_of_merged_list=True)
    del custom_gtf
    
    transcript_ids_wo_cds, NMD_features_df =\
    handle_cds_transcripts(transcript_gtftk_object, transcript_ids, NMD_features_df)
    

    start_no_cds = time.time() - start_time
    print("Time for start and checking no cds", start_no_cds)
    #get fasta of transcripts with known id
    if len(transcript_ids_wo_cds) > 0:
        print("known_tids_no_cds", len(transcript_ids_wo_cds))
        transcripts_with_CDS = determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, sys.argv[2], sys.argv[3], sys.argv[1])
        print(type(transcripts_with_CDS['end_ORF']))
        print(type(transcripts_with_CDS['last_exon_length']))
        transcripts_with_CDS.set_index('tid', inplace = True)
        transcripts_with_CDS['distance_stop_EJC'] = \
        transcripts_with_CDS['t_length'].astype(int) -\
        transcripts_with_CDS['last_exon_length'].astype(int) -\
        transcripts_with_CDS['end_ORF'].astype(int)

        transcripts_with_CDS['50_nt_rule'] = transcripts_with_CDS['distance_stop_EJC'].ge(50).astype(int)
        print(transcripts_with_CDS.head())
        NMD_features_df = NMD_features_df.join(transcripts_with_CDS, how='outer')
        NMD_features_df['50_nt'].fillna(NMD_features_df['50_nt_rule'])
        print(NMD_features_df.head())
    print(sum(NMD_features_df["50_nt_rule"].notna()), sum(NMD_features_df["50_nt_rule"][NMD_features_df["50_nt_rule"].notna()]))
    NMD_features_df.to_csv("NMD_features_df_NMD_trans.csv")
    return  0





if __name__ == "__main__":
    main()
