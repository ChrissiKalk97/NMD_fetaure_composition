import regex as re
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import itertools


def get_stop_codon_identity(CDS_seqs, NMD_features_df):
    """code the stop codon identity as boolean :'stop_TGA' or 
    'stop_TAA' or if none of the two it is TAG """

    NMD_features_df['stop_TGA'] = 0
    NMD_features_df['stop_TAA'] = 0
    for seq in CDS_seqs:
        if str(seq.seq[-3:]) == 'TGA':
            NMD_features_df.loc[seq.id.split(':')[0], 'stop_TGA'] = 1
        elif str(seq.seq[-3:]) == 'TAA':
            NMD_features_df.loc[seq.id.split(':')[0], 'stop_TAA'] = 1
    return NMD_features_df



def get_base_after_stop(transcript_sequences, NMD_features_df):
    """code the identity of the base after the stop codon """

    NMD_features_df['4th_stop_C'] = 0
    NMD_features_df['4th_stop_G'] = 0
    NMD_features_df['4th_stop_T'] = 0
    for seq in transcript_sequences:
        end_CDS = NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF']
        three_prime = str(seq.seq[int(end_CDS):])
        if len(three_prime) > 3:
            if three_prime[3] == 'C':
                NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_C'] = 1
            elif three_prime[3] == 'G':
                NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_G'] = 1
            elif three_prime[3] == 'T':
                NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_T'] = 1
        else:
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_C'] = pd.NA
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_G'] = pd.NA
            NMD_features_df.loc[seq.id.split(':')[0], '4th_stop_T'] = pd.NA
    return NMD_features_df



def get_GC_content_in30bp_ribo_window(transcript_sequences, NMD_features_df):
    """calculate the GC content in a window of 30bp centered around the stop codon"""

    NMD_features_df['GC_perc_30_bp_round_stop'] = 0.0
    for seq in transcript_sequences:
        end_CDS = int(NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF'])
        window_30 = str(seq.seq[end_CDS-15:end_CDS+15])
        C_count = window_30.count('C')
        G_count = window_30.count('G')
        NMD_features_df.loc[seq.id.split(':')[0], 'GC_perc_30_bp_round_stop'] = \
            (C_count + G_count) / 30
    return NMD_features_df



def get_GC_content_in15bp_ribo_window(transcript_sequences, NMD_features_df):
    """calculate the the GC content in both a 15bp window upstream and
    downstream of the stop codon"""

    NMD_features_df['GC_perc_up_15_bp_stop'] = 0.0
    NMD_features_df['GC_perc_down_15_bp_stop'] = 0.0
    for seq in transcript_sequences:
        end_CDS = int(NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF'])
        window_up_15 = str(seq.seq[end_CDS-15:end_CDS])
        window_down_15 = str(seq.seq[end_CDS:end_CDS+15])
        C_count_up = window_up_15.count('C')
        G_count_up = window_up_15.count('G')
        C_count_down = window_down_15.count('C')
        G_count_down = window_down_15.count('G')
        NMD_features_df.loc[seq.id.split(':')[0], 'GC_perc_up_15_bp_stop'] = \
            (C_count_up + G_count_up) / 15
        NMD_features_df.loc[seq.id.split(':')[0], 'GC_perc_down_15_bp_stop'] = \
            (C_count_down + G_count_down) / 15
    return NMD_features_df



def get_number_of_exons_transcript(transcript_sequences, NMD_features_df):
    """calculate number of exons in transcript and in the 3'UTR"""

    NMD_features_df['nr_exons_in_transcript'] = 0
    NMD_features_df['nr_exons_in_3prime'] = 0
    for seq in transcript_sequences:
        description = seq.description.split(':')
        NMD_features_df.loc[seq.id.split(':')[0],'nr_exons_in_transcript'] = \
        len(description) - 2
        length_3prime = NMD_features_df.loc[seq.id.split(':')[0], '3_UTR_length']
        exon_counter = 0
        if description[1] == 'strand-':
            while length_3prime > 0:
                exon_coords = description[2 + exon_counter].split('-')
                length_3prime = length_3prime - (int(exon_coords[2]) - int(exon_coords[1]) +1)
                exon_counter = exon_counter + 1
        else:
            while length_3prime > 0:
                #for plus strand get exons from the back
                exon_coords = description[len(description) - 1 - exon_counter].split('-')
                length_3prime = length_3prime - (int(exon_coords[2]) - int(exon_coords[1]) +1)
                exon_counter = exon_counter + 1
        NMD_features_df.loc[seq.id.split(':')[0],'nr_exons_in_3prime'] = \
        exon_counter
    return NMD_features_df



def find_UPF1_motifs_in3prime(transcript_sequences, NMD_features_df):
    """calculate absolute and relative occurrence of UPF1 motifs in the
    3'UTR"""

    NMD_features_df['UPF1_motifs_in3prime_total'] = 0
    NMD_features_df['UPF1_motifs_in3prime_relative'] = 0.0
    for seq in transcript_sequences:
        end_CDS = int(NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF'])
        three_prime = str(seq.seq[end_CDS:])
        UPF1_motif_count = three_prime.count('CTGGG')
        UPF1_motif_count = UPF1_motif_count + three_prime.count('CTGTG')
        NMD_features_df.loc[seq.id.split(':')[0],'UPF1_motifs_in3prime_total'] =\
        UPF1_motif_count
        NMD_features_df.loc[seq.id.split(':')[0],'UPF1_motifs_in3prime_relative'] =\
        UPF1_motif_count/len(three_prime)
    return NMD_features_df




def count_k_mers(transcript_sequences, NMD_features_df):
    """Count the amount of DNA k-mers in a window of 30bp
      centered around the stop codon"""
    
    # Define the DNA bases
    dna_bases = ['A', 'C', 'G', 'T']
    # Generate all possible 4-mers
    all_4mers = [''.join(x) for x in itertools.product(dna_bases, repeat=4)]
    for k_mer in all_4mers:
        NMD_features_df[k_mer] = 0
        for seq in transcript_sequences:
            end_CDS = int(NMD_features_df.loc[seq.id.split(':')[0], 'end_ORF'])
            window_30 = str(seq.seq[end_CDS-15:end_CDS+15])
            NMD_features_df.loc[seq.id.split(':')[0], k_mer]\
                  = window_30.count(k_mer)
    return NMD_features_df
