import time
import regex as re
from pybedtools import BedTool
import pybedtools as pb
import pandas as pd

from new_helper_functions import get_transcript_string
from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, get_ORF_start_by_gene, coinciding_start_sites, filter_bed_file
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds[0:1000]}
    genome_file = reference_sequences_file
    #"./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa"
    start_time = time.time()
    sequences = get_fasta_tid(transcripts_no_cds, genome_file)
    print("Number of transcripts for which Sequence was generated", len(set([rec.id for rec in sequences])))
    get_fasta = time.time() - start_time
    ORFs = OrfFinder(sequences)
    getORFs = time.time() - start_time - get_fasta 
    cds_bed_positions, start_positions = get_cds_genomic_coordinates(ORFs)
    genomic_time = time.time() - start_time - get_fasta - getORFs 
    

    #extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
    
    reference_genes, reference_gtf = get_pct_reference(reference_file, gene_ids_ORF_transcripts)
    print("nr of genes for which gene not found or no protein coding transcripts",\
           len([gene_id for gene_id in gene_ids_ORF_transcripts if gene_id not in reference_genes.keys()]))

    orf_start_sites_by_gene = get_ORF_start_by_gene(start_positions)
    nr_trans = {k:[t[4] for t in v] for (k,v) in orf_start_sites_by_gene.items()}
    nr_trans = [len(set(v)) for v in nr_trans.values()]
    print("nr of transcripts with calculated start positions", sum(nr_trans))


    transcripts_wo_pc_reference = {k:[t[4] for t in v] for (k,v) in orf_start_sites_by_gene.items() if k not in reference_genes.keys()}
    nr_tr_no_pc_ref = [len(set(v)) for v in transcripts_wo_pc_reference.values()]
    prepare_reference = time.time() - start_time - get_fasta - getORFs - genomic_time

    transcripts_cds_determined, transcripts_several_orfs =\
        coinciding_start_sites(gene_ids_ORF_transcripts, reference_genes, orf_start_sites_by_gene)
    t_co_ss = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference

    


    transcripts_several_orfs_bed = filter_bed_file(transcripts_several_orfs, cds_bed_positions)
    #this filtering step worked


    transcripts_several_orfs_string = get_transcript_string(transcripts_several_orfs)
    reference_bed = BedTool(reference_gtf.select_by_key("feature", "CDS")\
        .select_by_key("transcript_id", transcripts_several_orfs_string)\
            .to_bed(name=('gene_id', 'transcript_id'), sep=":"))
    
    protein_coding_reference_transcripts_bed = filter_bed_file(transcripts_several_orfs, reference_bed)
    t_filter_beds = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference - t_co_ss

    intersection = transcripts_several_orfs_bed.intersect(protein_coding_reference_transcripts_bed, wo = True).saveas()
    #sum intersection by name of file a: splitORFs
    #summed_counts = intersection.groupby(g=[4, 10], c=13, o="sum", output = "ORF_pct_intersection_summed.csv")
    #summed_counts = pd.read_csv("ORF_pct_intersection_summed.csv", sep ="\t")
    summed_counts = pd.read_table(intersection.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand',\
                                               'chrom_tar', 'start_tar', 'stop_tar', 'name_tar', 'score_tar', 'strand_tar', "overlap"])

    summed_counts = summed_counts.groupby(["name", "name_tar"])['overlap'].sum()
    summed_counts = summed_counts.reset_index()
    time_intersect_sum = time.time() - start_time - get_fasta - getORFs -\
          genomic_time - prepare_reference - t_co_ss - t_filter_beds
    


    #print(type(summed_counts["name"].str.split(":")))
    #print(summed_counts["name"].str.split(":"))
    #print(summed_counts["name"].str.split(":")[0])
    #extract transcript id from name and search target with largest overlap
    #the corresponding ORF will be the CDS for the transcript

    print("nr of transcripts for which gene not found or no protein coding reference", sum(nr_tr_no_pc_ref))
    print("Nr transcripts with several ORFs", len(transcripts_several_orfs))
    print("Nr transcripts with only one CDS", len(transcripts_cds_determined.keys()))

    pb.cleanup(remove_all=True)


    print("time needed for fasta", get_fasta)
    print("time needed for ORFs", getORFs)
    print("time needed for genomic", genomic_time)
    print("time needed for reference preparation", prepare_reference)
    print("time needed for coinc start sites", t_co_ss)
    print("Time needed to filter beds", t_filter_beds)
    print("Time needed to intersect and sum counts", time_intersect_sum)


    

    