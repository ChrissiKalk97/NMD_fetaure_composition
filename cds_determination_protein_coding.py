import pybedtools as pb
from pygtftk.gtf_interface import GTF
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates, \
    get_pct_reference, find_cds_orf, get_transcript_string
from OrfFinder_py3 import OrfFinder


def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds,
                  reference_file, reference_sequences_file, output_name):
    """determine CDS of transcripts in gtf file by checking for protein sequence
    overlap with known CDS from the same genomic region, then selecting the most 5'
    ORF"""

    # obtain transcript information for transcripts of which CDS needs to be determined
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in 
            transcript_gtftk_object.keys() if k in set(transcript_ids_wo_cds)}

    # define genome file
    genome_file = reference_sequences_file

    # get fasta sequence of input transcripts
    sequences = get_fasta_tid(transcripts_no_cds, genome_file, seq_type='exon')
    SeqIO.write(
        sequences, f'Output/{output_name}_sequences_transcripts_provided.fasta', 'fasta')

    # find all possible ORFs
    ORFs = OrfFinder(sequences)
    ORFs_for_fasta = [SeqRecord(id=str(orf.name) + '|' + str(orf.id),
                                seq=orf.seq.translate(), description='') for orf in ORFs]
    SeqIO.write(ORFs_for_fasta,
                f'Output/{output_name}_ORFFinder.fasta', 'fasta')
    del ORFs_for_fasta

    # extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))

    reference_gtf_CDS = get_pct_reference(
        reference_file, gene_ids_ORF_transcripts)

    # get possible targets from reference for the CDS overlap
    gene_string = get_transcript_string(gene_ids_ORF_transcripts)
    reference_tragets = GTF(reference_file, check_ensembl_format=False)\
        .select_by_key('gene_biotype', 'protein_coding')\
        .select_by_key('gene_id', gene_string)\
        .extract_data('transcript_id,start,end,exon_number,feature,\
                      strand,chrom,gene_id,score',
                      as_dict_of_merged_list=True)

    # get protein coding genes that are present in the custom gtf file
    protein_coding_genes = [target_info[6]
                            for target_info in reference_tragets.values()]
    protein_coding_genes = set(protein_coding_genes)

    # ORF seqeunces in DNA
    ORFs_for_fasta = [SeqRecord(id=str(orf.name) + '|' + str(orf.id),
                                seq=orf.seq, description='') for orf in ORFs]
    SeqIO.write(ORFs_for_fasta,
                f'Output/{output_name}_ORFFinder_DNA.fasta', 'fasta')
    del ORFs_for_fasta

    # filter for protein coding genes
    ORFs = [ORF for ORF in ORFs if ORF.name in protein_coding_genes]
    SeqIO.write(
        ORFs, f'Output/{output_name}_ORFS_protein_coding_genes.fasta', 'fasta')

    # get cds coordinates genomic
    orf_bed_positions, orf_dict_exon_with_stop_length, distance_stop_next_EJC\
        = get_cds_genomic_coordinates(ORFs)
    del ORFs

    # obtain sequences of the target transcripts for exact matching
    target_sequences = get_fasta_tid(
        reference_tragets, genome_file, seq_type='CDS')
    target_sequences = [SeqRecord(id=str(target.name) + '|' + str(target.id),
                                  seq=target.seq.translate(), description='') 
                                  for target in target_sequences]
    SeqIO.write(target_sequences,
                f'Output/{output_name}_target_sequences.fasta', 'fasta')
    del target_sequences

    # determine the CDS for the transcripts
    transcripts_with_CDS = find_cds_orf(reference_gtf_CDS, 
                                        orf_bed_positions,
                                        f'Output/{output_name}_ORFFinder.fasta', 
                                        f'Output/{output_name}_target_sequences.fasta')

    # clean-up
    pb.cleanup(remove_all=True)

    # filter for transcripts, for which an ORFs was found
    transcripts_of_interest = transcripts_with_CDS['name'].str.split(
        r':|\|').str[1]
    sequences = [
        seq for seq in sequences if seq.id in transcripts_of_interest.values]

    SeqIO.write(
        sequences, f'Output/{output_name}_transcripts_of_interest.fasta', 'fasta')

    # map the exon length of the exon with stop to the table
    transcripts_with_CDS['exon_with_stop_length'] =\
        transcripts_with_CDS['name'].map(orf_dict_exon_with_stop_length)

    # map the exon length of the exon with stop to the table
    transcripts_with_CDS['distance_stop_next_EJC'] =\
        transcripts_with_CDS['name'].map(distance_stop_next_EJC)

    return transcripts_with_CDS, sequences
