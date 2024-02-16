import regex as re
from pybedtools import BedTool

from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, get_ORF_start_by_gene, coinciding_start_sites, filter_bed_file
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds[0:10000]}
    genome_file = reference_sequences_file
    #"./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa"
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
    transcripts_wo_pc_reference = {k:[t[4] for t in v] for (k,v) in orf_start_sites_by_gene.items() if k not in reference_genes.keys()}
    nr_tr_no_pc_ref = [len(set(v)) for v in transcripts_wo_pc_reference.values()]
    

    transcripts_cds_determined, transcripts_several_orfs =\
        coinciding_start_sites(gene_ids_ORF_transcripts, reference_genes, orf_start_sites_by_gene)
    
    transcripts_several_orfs_bed = filter_bed_file(transcripts_several_orfs, cds_bed_positions)
    #this filtering step worked

    reference_bed = BedTool(reference_file)
    protein_coding_reference_transcripts_bed = filter_bed_file(transcripts_several_orfs, reference_bed)
    print(protein_coding_reference_transcripts_bed)
    print(type(protein_coding_reference_transcripts_bed))
    intersection = transcripts_several_orfs_bed.intersect(protein_coding_reference_transcripts_bed, wa = True, f = 1.0).saveas()
    #print()

    print("nr of ORFs for which gene not found or no protein coding reference", sum(nr_tr_no_pc_ref))
    print("Nr transcripts with several ORFs", len(transcripts_several_orfs))
    print("Nr transcripts with only one CDS", len(transcripts_cds_determined.keys()))

    #pybedtools.cleanup(remove_all=True)

    #filter ORFs from ORF Finder for CDS that are determined, return