#######################################################################
##undefined_CDS_transcripts.py counts the number of transcripts in a ##
##gtf file which have any of the supplied tags                       ##
##then it checks if a csv fiel output of the fifty nt rule is provided
##how many of the entries with no info belong to this group         ###
#######################################################################
##usage

#python undefined_CDS_transcripts.py gtf.gtf tag1,tag2,tag3 out.csv



import sys
from typing import List
import pandas as pd
from pygtftk.gtf_interface import GTF

def get_transcript_string(transcript_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    transcript_string = ''
    if len(transcript_ids) > 1: 
        for transcript_id in transcript_ids[:-1]:
            transcript_string += transcript_id+','
    transcript_string += transcript_ids[-1]
    return transcript_string



#GTF file to check
custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

#get all tags
tags = sys.argv[2]#'mRNA_start_NF,cds_end_NF,cds_start_NF,mRNA_end_NF'

#read the output of the 50nt rule
fifty_out = pd.read_csv(sys.argv[3], index_col=0)
print(fifty_out.head())


custom_gtf_NF = custom_gtf.select_by_key('feature', 'transcript')\
.select_by_key('tag', tags)
unique_transcripts_NF = custom_gtf_NF.get_tx_ids(nr=True)
print('number of transcripts with undefined CDS/mRNA end/start: ', len(unique_transcripts_NF))

#filter the gtf file such that the NF entries are filtered out
all_tids = custom_gtf.get_tx_ids(nr=True)
filtered_tids = [tid for tid in all_tids if tid not in unique_transcripts_NF]
filtered_tids_string = get_transcript_string(filtered_tids)
custom_gtf.select_by_key('feature', 'transcript,exon')\
.select_by_key('transcript_id', filtered_tids_string).write(sys.argv[5])


not_found_fifty = fifty_out[fifty_out['name_tar'].isnull()].index
#print(fifty_out[fifty_out['name_tar'].isnull()].head())
print('number of transcripts for which no target found with 50nt rule: ', len(not_found_fifty))

not_found_and_NF = [trans for trans in not_found_fifty if trans in unique_transcripts_NF]
print('number of transcripts for which no target found with 50nt rule and which are NF: ', len(not_found_and_NF))

#print([trans for trans in not_found_fifty if trans not in unique_transcripts_NF])
wrong = float(sys.argv[4])
fifty_out['50_nt'] = fifty_out['50_nt'].astype(float)
wrong_class = fifty_out[fifty_out['50_nt'] == wrong]
wrong_class_tids = wrong_class.index
wrong_class_and_nf = [trans for trans in wrong_class_tids if trans in unique_transcripts_NF]

print("number of transcripts wrongly classified: ", len(wrong_class_tids))
print("number of transcripts wrongly classified and NF: ", len(wrong_class_and_nf))

