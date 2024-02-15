import regex as re

from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, get_ORF_start_by_gene, coinciding_start_sites
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds[0:100000]}
    genome_file = reference_sequences_file#"./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa"
    sequences = get_fasta_tid(transcripts_no_cds, genome_file)
    ORFs = OrfFinder(sequences)
    cds_bed_positions, start_positions = get_cds_genomic_coordinates(ORFs)
    
    #extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
    
    reference_genes = get_pct_reference(reference_file, gene_ids_ORF_transcripts)
    print("gene not found or no protein coding transcripts",\
           len([gene_id for gene_id in gene_ids_ORF_transcripts if gene_id not in reference_genes.keys()]))

    orf_start_sites_by_gene = get_ORF_start_by_gene(start_positions)

    transcripts_cds_determined, bed_for_intersection =\
        coinciding_start_sites(gene_ids_ORF_transcripts, reference_genes, orf_start_sites_by_gene, cds_bed_positions)
    
    print(bed_for_intersection)