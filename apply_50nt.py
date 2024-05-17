from operator import itemgetter
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from helper_functions_cds_determination import get_transcript_string

def get_length_last_exon(transcript_ids_list, gtf_file):
    '''extract last exon length per transcript from reference'''
    last_exon_length_by_transcript = {}
    transcript_string = get_transcript_string(transcript_ids_list)
    transcript_entries = gtf_file\
        .select_by_key('transcript_id', transcript_string)\
        .select_by_key('feature', 'exon')\
        .extract_data('transcript_id,start,end,exon_number,feature,strand',
                       as_dict_of_merged_list=True)
    for transcript_id in transcript_ids_list:
        transcript_exons = transcript_entries[transcript_id]
        transcript_exons = [transcript_exons[x:x+5] for x in range(0, len(transcript_exons), 5)]
        strand = transcript_exons[0][4]
        #sort transcript exons in ascending order
        transcript_exons = sorted(transcript_exons, key = itemgetter(int(0)))
        if strand == '+':
            last_exon = transcript_exons[-1]
        else:
            last_exon = transcript_exons[0]
        last_exon_length_by_transcript[transcript_id] = int(last_exon[1]) - int(last_exon[0]) + 1
    return last_exon_length_by_transcript


def calculate_50nt_rule(transcripts_with_CDS: pd.DataFrame, sequences) -> pd.DataFrame:
    t_length_dict = {}
    for sequence_rec in sequences:
        t_length_dict[sequence_rec.id] = len(sequence_rec.seq)
    
    transcripts_with_CDS['t_length'] = transcripts_with_CDS['tid'].map(t_length_dict)
    transcripts_with_CDS.set_index('tid', inplace = True)
    transcripts_with_CDS['distance_stop_EJC'] = \
    transcripts_with_CDS['t_length'].astype(int) -\
    transcripts_with_CDS['last_exon_length'].astype(int) -\
    transcripts_with_CDS['end_ORF']
    transcripts_with_CDS['50_nt_rule'] = transcripts_with_CDS['distance_stop_EJC'].ge(50).astype(int)
    return transcripts_with_CDS
