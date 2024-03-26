import time
import regex as re
from pybedtools import BedTool
import pybedtools as pb
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, find_cds_orf, get_length_last_exon,\
      calculate_50nt_rule, get_transcript_string
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file, gtf_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds}
    genome_file = reference_sequences_file
    #'./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa'
    start_time = time.time()
    sequences = get_fasta_tid(transcripts_no_cds, genome_file, seq_type = 'exon')
    SeqIO.write(sequences, f'Output/sequences_transcripts_provided_integer.fasta', 'fasta')
    get_fasta = time.time() - start_time
    ORFs = OrfFinder(sequences)
    ORFs_for_fasta = [SeqRecord(id = str(orf.name) + '|' + str(orf.id),\
                                   seq = orf.seq.translate(), description = '') for orf in ORFs]
    SeqIO.write(ORFs_for_fasta, f'Output/ORFFinder_integer.fasta', 'fasta')

    #transcript sequences w/o ORFs:
    transcripts_with_ORFs = [trans.id for trans in ORFs_for_fasta]
    transcripts_no_orfs = [trans for trans in sequences if trans.id not in transcripts_with_ORFs]
    SeqIO.write(sequences, f'Output/transcripts_w_o_ORFs_Ensmebl.fasta', 'fasta')


    ORFs_for_fasta = [SeqRecord(id = str(orf.name) + '|' + str(orf.id),\
                                   seq = orf.seq, description = '') for orf in ORFs]
    SeqIO.write(ORFs_for_fasta, f'Output/ORFFinder_DNA.fasta', 'fasta') 



    getORFs = time.time() - start_time - get_fasta 
    orf_bed_positions = get_cds_genomic_coordinates(ORFs)
    genomic_time = time.time() - start_time - get_fasta - getORFs 
    

    #extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
    
    reference_gtf_CDS = get_pct_reference(reference_file, gene_ids_ORF_transcripts)
    prepare_reference = time.time() - start_time - get_fasta - getORFs - genomic_time

    #get target sequences as a fasta file
    gene_string = get_transcript_string(gene_ids_ORF_transcripts)
    target_transcript_ids = GTF(reference_file, check_ensembl_format=False)\
        .select_by_key('feature', 'CDS')\
        .select_by_key('gene_biotype', 'protein_coding')\
        .select_by_key('gene_id', gene_string)\
        .extract_data('transcript_id,start,end,exon_number,feature,strand,chrom,gene_id,score',
                       as_dict_of_merged_list=True)
    
    target_sequences = get_fasta_tid(target_transcript_ids, genome_file, seq_type = 'CDS')
    target_sequences = [SeqRecord(id = str(target.name) + '|' + str(target.id),\
                                   seq = target.seq.translate(), description = '') for target in target_sequences]
    print('target_sequences', target_sequences[:5])
    SeqIO.write(target_sequences, f'Output/target_sequences_integer.fasta', 'fasta')
    #reference_gtf = GTF(reference_file, check_ensembl_format=False)
    #gene_string = get_transcript_string(gene_ids_ORF_transcripts)    
    #reference_gtf_CDS = reference_gtf\
    #.select_by_key('gene_id', gene_string)\
    #.select_by_key('gene_biotype', 'protein_coding')
    
    transcripts_with_CDS  = find_cds_orf(reference_gtf_CDS, orf_bed_positions, 'Output/ORFFinder_integer.fasta', 'Output/target_sequences_integer.fasta')
    t_find_cds_orf = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference

    print('transcripts with their partners', transcripts_with_CDS.head())
    last_exon_length_dict = get_length_last_exon(transcripts_with_CDS['tid'].to_list(), gtf_file)
    #print(len(last_exon_length_dict) == len(transcripts_with_CDS.index))
    transcripts_with_CDS['last_exon_length'] = transcripts_with_CDS['tid'].map(last_exon_length_dict)
    

    transcripts_with_CDS = calculate_50nt_rule(transcripts_with_CDS, sequences)
    #get the transcript length
    
    pb.cleanup(remove_all=True)


    print('time needed for fasta', get_fasta)
    print('time needed for ORFs', getORFs)
    print('time needed for genomic', genomic_time)
    print('time needed for reference preparation', prepare_reference)
    print('time needed for coinc finding ORF for CDS', t_find_cds_orf)
    


    return transcripts_with_CDS

    