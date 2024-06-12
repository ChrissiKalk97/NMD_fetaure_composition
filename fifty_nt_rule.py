# fifty_nt_rule.py determines the CDS for given transcripts in a gtf
# (if not supplied), applies the 50 nt rule and calculates other sequence
# features of the transcripts
# the output is a csv file with the computed feature set for each transcript
# with information of the chosen CDS (in case this was not present in the gtf file)

# uasge: python fifty_nt_rule.py custom.gtf Ensembl.gtf genome.fasta output_name.csv

import os
import sys
import time

import numpy as np
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
from warnings import simplefilter


from CDS_annotation_present import handle_cds_transcripts
from cds_determination_protein_coding import determine_cds
from apply_50nt import get_length_last_exon, calculate_50nt_rule
from helper_functions_cds_determination import get_fasta_tid
from calculate_sequence_features import get_stop_codon_identity, get_base_after_stop, \
    get_GC_content_in30bp_ribo_window, get_number_of_exons_transcript, find_UPF1_motifs_in3prime, \
    get_GC_content_in15bp_ribo_window, count_k_mers, optimal_codon_usage


def main():
    
    # create Output directory
    output_name = sys.argv[4].split('/')[-1][:-4]
    folder_path = f'Output/{output_name}'
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
        print(f"Folder '{folder_path}' created.")
    else:
        print(f"Folder '{folder_path}' already exists.")
    # create correct output_name
    output_name = f'{output_name}/{output_name}'

    # ignore pandas warnings
    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    #################################################################################
    ### start with the feature calculation for the transcriptome assemblies###########
    #################################################################################

    # read in gtf file
    custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

    # get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    transcript_ids = sorted(transcript_ids)
    print('number of transcripts in custom gtf: ', len(transcript_ids))

    # build dataframe to store the computed features
    NMD_features_df = pd.DataFrame(columns=['50_nt', 'last_exon_length',
                                            't_length', 'start_ORF', 'end_ORF',
                                            'exon_with_stop_length',
                                            'stop_TGA', 'stop_TAA'],
                                   index=transcript_ids)

    # get exons, CDS and stop codons from the gtf object as dict
    transcript_gtftk_object = custom_gtf\
        .select_by_key('feature', 'CDS,exon,stop_codon')\
        .extract_data('transcript_id,start,end,exon_number,feature,strand,chrom,gene_id,score',
                      as_dict_of_merged_list=True)

    #################################################################################
    ### handling of transcripts with CDS anno in custom gtf###########################
    ### note down transcript ids for which no CDS anno found##########################
    #################################################################################

    transcript_ids_wo_cds, NMD_features_df =\
        handle_cds_transcripts(transcript_gtftk_object,
                               transcript_ids, NMD_features_df)

    # get transcript and ORF sequences for the transcripts with CDS
    transcript_ids_wo_cds_set = set(transcript_ids_wo_cds)
    tids_with_cds = [
        tid for tid in transcript_ids if tid not in transcript_ids_wo_cds_set]
    tids_with_cds_set = set(tids_with_cds)
    CDS_seqs = []
    transcript_seqs = []
    if len(tids_with_cds_set) > 0:
        transcripts_with_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()
                                if k in tids_with_cds_set}
        transcript_seqs = get_fasta_tid(
            transcripts_with_cds, sys.argv[3], seq_type='exon')
        CDS_seqs = get_fasta_tid(transcripts_with_cds,
                                 sys.argv[3], seq_type='CDS', plus_stop=True)

    #################################################################################
    ### handle transcripts with no CDS annoation: CDS determintation, then 50 nt rule#
    #################################################################################
    if len(transcript_ids_wo_cds) > 0:
        print('transcripts for which no CDS annotation was given in the custom gtf: ',
              len(transcript_ids_wo_cds))
        time_0 = time.time()
        # determine CDS for source transcripts
        transcripts_calculated_CDS, sequences = determine_cds(transcript_gtftk_object, transcript_ids_wo_cds,
                                                              sys.argv[2], sys.argv[3], output_name)
        CDS_time = time.time() - time_0
        print("CDS_time_required:", CDS_time, flush=True)

        # drop all for which no ORF was determined
        transcripts_calculated_CDS = transcripts_calculated_CDS.dropna(
            subset='name')

        time_1 = time.time()
        # get the length of the last exon per transcript for 50 nt rule (from the custom gtf)
        last_exon_length_dict = get_length_last_exon(
            transcripts_calculated_CDS['tid'].to_list(), custom_gtf)
        transcripts_calculated_CDS['last_exon_length'] = transcripts_calculated_CDS['tid'].map(
            last_exon_length_dict)

        # apply 50nt rule
        transcripts_calculated_CDS = calculate_50nt_rule(
            transcripts_calculated_CDS, sequences)
        
        fifty_nt_time = time.time()-time_1

        # select the sequences of the ORFs selected as the CDS
        transcripts_calculated_CDS['ORF_id'] = transcripts_calculated_CDS['name'].str.split(
            '|').str[1]
        transcripts_calculated_CDS['ORF_id'] = transcripts_calculated_CDS['ORF_id'].astype(
            'string')
        ORFs = SeqIO.parse(
            f'Output/{output_name}_ORFS_protein_coding_genes.fasta', "fasta")
        ORF_list = []
        ORF_id_set = set(transcripts_calculated_CDS['ORF_id'].unique())
        for ORF in ORFs:
            if ORF.id in ORF_id_set:
                ORF_list.append(ORF)
        ORFs = ORF_list

        prepare_ORFs_time = time.time()-time_1-fifty_nt_time

        # combine CDS sequences for given and calculated CDS
        CDS_seqs = CDS_seqs + ORFs
        CDS_time= time.time()-time_1-fifty_nt_time-prepare_ORFs_time 

        # join df of calculated CDS transcripts with df of the CDS annotated transcripts
        NMD_features_df = pd.concat(
            [transcripts_calculated_CDS, NMD_features_df], axis=0, join='outer')
        NMD_features_df['50_nt'] = np.where(NMD_features_df['50_nt'].isna(),
                                            NMD_features_df['50_nt_rule'], NMD_features_df['50_nt'])
        NMD_features_df.drop(['50_nt_rule'], axis=1, inplace=True)

        # combine transcript sequences with given and calculated CDS
        transcript_seqs = transcript_seqs + sequences
        concat_time = time.time()-time_1-fifty_nt_time-prepare_ORFs_time-CDS_time

    # drop rows with NA's, for those transcripts no CDS could be determined
    NMD_features_df = NMD_features_df.dropna(subset=['t_length', 'start_ORF'])

    # calculate features used for the NMD classifier
    NMD_features_df = get_stop_codon_identity(CDS_seqs, NMD_features_df)
    NMD_features_df['3_UTR_length'] = NMD_features_df['t_length'] - \
        NMD_features_df['end_ORF']
    NMD_features_df['5_UTR_length'] = NMD_features_df['start_ORF']
    NMD_features_df['distance_stop_from_start'] = NMD_features_df['end_ORF'] - \
        NMD_features_df['start_ORF'] - 2
    NMD_features_df['stop_150bp_from_start'] = np.where(
        NMD_features_df['distance_stop_from_start'] > 150, 0, 1)
    NMD_features_df = get_base_after_stop(transcript_seqs, NMD_features_df)
    NMD_features_df = get_GC_content_in30bp_ribo_window(
        transcript_seqs, NMD_features_df)
    NMD_features_df = get_GC_content_in15bp_ribo_window(
        transcript_seqs, NMD_features_df)
    NMD_features_df = get_number_of_exons_transcript(
        transcript_seqs, NMD_features_df)
    NMD_features_df = find_UPF1_motifs_in3prime(
        transcript_seqs, NMD_features_df)
    NMD_features_df = count_k_mers(transcript_seqs, NMD_features_df)
    NMD_features_df = optimal_codon_usage(transcript_seqs, NMD_features_df)
    
    feature_time = time.time()-time_1-fifty_nt_time-prepare_ORFs_time-CDS_time-concat_time
    print("50nt_time_required:", fifty_nt_time, flush=True)
    print("prepare_ORFs_time_required:", prepare_ORFs_time, flush=True)
    print("CDS_prep_time_required:", CDS_time, flush=True)
    print("concat_time_required:", concat_time, flush=True)
    print("feature_time_required:", feature_time, flush=True)

    # print output
    print(NMD_features_df.head())
    print('number of transcripts for which 50 nt rule was calculated: ',
          sum(NMD_features_df['50_nt'].notna()))
    print('number of transcripts found to be NMD sensitive (no escape) by 50nt rule ',
          sum(NMD_features_df['50_nt'][NMD_features_df['50_nt'].notna()]))

    # drop unnecessary columns
    if len(transcript_ids_wo_cds) > 0:
        NMD_features_df.drop(['gid', 'gid_target', 'ORF_nr', 'target_length', 'name',
                              'target_coverage_percentage', 'protein_overlap_perc', 'overlap',
                              'protein_overlap_aa'], axis=1, inplace=True)

    # write results to csv
    NMD_features_df.to_csv(f'Output/{output_name}.csv')
    return 0


if __name__ == '__main__':
    main()
