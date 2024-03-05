import time
import regex as re
from pybedtools import BedTool
import pybedtools as pb
import pandas as pd

from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, find_cds_orf, get_length_last_exon,\
      calculate_50nt_rule
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file, gtf_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds}
    genome_file = reference_sequences_file
    #'./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa'
    start_time = time.time()
    sequences = get_fasta_tid(transcripts_no_cds, genome_file)
    get_fasta = time.time() - start_time
    ORFs = OrfFinder(sequences)
    #for orf in ORFs:
        #print(orf.id, orf.seq)
    getORFs = time.time() - start_time - get_fasta 
    orf_bed_positions = get_cds_genomic_coordinates(ORFs)
    genomic_time = time.time() - start_time - get_fasta - getORFs 
    

    #extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
    
    reference_gtf_CDS = get_pct_reference(reference_file, gene_ids_ORF_transcripts)
    prepare_reference = time.time() - start_time - get_fasta - getORFs - genomic_time

    
    transcripts_with_CDS  = find_cds_orf(reference_gtf_CDS, orf_bed_positions)
    t_find_cds_orf = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference

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

    