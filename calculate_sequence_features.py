import regex as re
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_stop_codon_identity(CDS_seqs, NMD_features_df):
    NMD_features_df['stop_TGA'] = 0
    NMD_features_df['stop_TAA'] = 0
    for seq in CDS_seqs:
        if str(seq.seq[-3:]) == 'TGA':
            NMD_features_df.loc[seq.id.split(':')[0], 'stop_TGA'] = 1
        elif str(seq.seq[-3:]) == 'TAA':
            NMD_features_df.loc[seq.id.split(':')[0], 'stop_TAA'] = 1
    return NMD_features_df

def get_base_after_stop(transcript_sequences, NMD_features_df):
    #NMD_features_df['start_ORF'] and NMD_features_df['end_ORF']
    #will give the CDS coordinates
    NMD_features_df['4th_stop_C'] = 0
    NMD_features_df['4th_stop_G'] = 0
    NMD_features_df['4th_stop_T'] = 0
    for seq in transcript_sequences:
        end_CDS = NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF']
        three_prime = str(seq.seq[end_CDS:])
        print(seq.description)
        if three_prime[4] == 'C':
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_C'] = 1
        elif three_prime[4] == 'G':
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_G'] = 1
        elif three_prime[4] == 'T':
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_T'] = 1
    return NMD_features_df

def get_GC_content_in30bp_ribo_window(transcript_sequences, NMD_features_df):
    NMD_features_df['GC_perc_30_bp_round_stop'] = 0.0
    for seq in transcript_sequences:
        end_CDS = NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF']
        window_30 = str(seq.seq[end_CDS-15:end_CDS+15])
        C_count = window_30.count('C')
        G_count = window_30.count('G')
        NMD_features_df.loc[seq.id.split(':')[0], 'GC_perc_30_bp_round_stop'] = \
            (C_count + G_count) / 30
    return NMD_features_df
