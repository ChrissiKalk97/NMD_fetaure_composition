import time
import regex as re
from pybedtools import BedTool
import pybedtools as pb
import pandas as pd

from new_helper_functions import get_transcript_string
from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, get_ORF_start_by_gene, coinciding_start_sites, filter_bed_file, get_length_last_exon
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
    getORFs = time.time() - start_time - get_fasta 
    cds_bed_positions, start_positions = get_cds_genomic_coordinates(ORFs)
    genomic_time = time.time() - start_time - get_fasta - getORFs 
    

    #extract gene ids of the transcripts that make ORFs
    gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
    gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
    
    reference_genes, reference_gtf = get_pct_reference(reference_file, gene_ids_ORF_transcripts)
    #print('nr of genes for which gene not found or no protein coding transcripts',\
      #     len([gene_id for gene_id in gene_ids_ORF_transcripts if gene_id not in reference_genes.keys()]))

    orf_start_sites_by_gene = get_ORF_start_by_gene(start_positions)
    nr_trans = {k:[t[4] for t in v] for (k,v) in orf_start_sites_by_gene.items()}
    nr_trans = [len(set(v)) for v in nr_trans.values()]
    #print('nr of transcripts with calculated start positions', sum(nr_trans))


    transcripts_wo_pc_reference = {k:[t[4] for t in v] for (k,v)\
            in orf_start_sites_by_gene.items() if k not in reference_genes.keys()}
    nr_tr_no_pc_ref = [len(set(v)) for v in transcripts_wo_pc_reference.values()]
    prepare_reference = time.time() - start_time - get_fasta - getORFs - genomic_time

    transcripts_cds_determined, transcripts_several_orfs, genes_tids_several_orfs =\
        coinciding_start_sites(gene_ids_ORF_transcripts, reference_genes, orf_start_sites_by_gene)
    t_co_ss = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference

    print('nr of transcripts for which gene not found or no protein coding reference', sum(nr_tr_no_pc_ref))
    print('Nr transcripts with several ORFs', len(transcripts_several_orfs))
    print('Nr transcripts with only one CDS', len(transcripts_cds_determined))


    #filter cds genomic positions of calculated ORFs for the transcripts that
    #have several ORFs to consider as CDS
    transcripts_several_orfs_bed = filter_bed_file(transcripts_several_orfs, cds_bed_positions).saveas('transcripts_with_several_orfs.bed')
    
    #get all transcript ids for which several ORFs are considered
    gene_several_orfs_string = get_transcript_string(list(set(genes_tids_several_orfs)))
    #create reference bed file containing CDS sites of protein coding transcripts
    #of the transcripts for which the CDS needs to be determined
    pc_reference_bed = BedTool(reference_gtf.select_by_key('feature', 'CDS')\
        .select_by_key('gene_id', gene_several_orfs_string)\
        .select_by_key('transcript_biotype', 'protein_coding')\
            .to_bed(name=('gene_id', 'transcript_id'), sep=':')).saveas("pc_reference.bed")
    #pb.cleanup(remove_all=True)


    t_filter_beds = time.time() - start_time - get_fasta - getORFs - genomic_time - prepare_reference - t_co_ss
    intersection_v = transcripts_several_orfs_bed.intersect(pc_reference_bed, v = True).saveas('intersection_v.bed')
    intersection_v = pd.read_table(intersection_v.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    intersection_v['tid'] = intersection_v['name'].str.extract(r'[A-Z0-9\.]*:([A-Z0-9\.]*):.*')


    intersection = transcripts_several_orfs_bed.intersect(pc_reference_bed, wao = True).saveas('intersection.bed')
    summed_counts = pd.read_table(intersection.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand',\
                    'chrom_tar', 'start_tar', 'stop_tar', 'name_tar', 'score_tar', 'strand_tar', 'overlap'])
    
    summed_counts = summed_counts.groupby(['name', 'name_tar'])['overlap'].sum()
    summed_counts = summed_counts.reset_index()
    time_intersect_sum = time.time() - start_time - get_fasta - getORFs -\
          genomic_time - prepare_reference - t_co_ss - t_filter_beds
    
    summed_counts['tid'] = summed_counts['name'].str.extract(r'[A-Z0-9\.]*:([A-Z0-9\.]*):.*')
    rowIds = summed_counts.groupby('tid')['overlap'].idxmax()
    transcripts_with_CDS = summed_counts.loc[rowIds]
    
    #calculate tids of protein coding genes, but without overlap with protein coding transcripts
    tid_no_overlap =\
        [tid for tid in intersection_v['tid'].to_list() if tid not in transcripts_with_CDS['tid'].to_list()]
    print("tid_no_overlap", set(tid_no_overlap))

    transcripts_with_CDS = transcripts_with_CDS.reset_index()[["name", "name_tar", "tid"]]
    #append the others with only one CDS 
    transcripts_cds_determined = pd.DataFrame(transcripts_cds_determined)
    
    print("length only one cds", len(transcripts_cds_determined.index))
    print("length several cds", len(transcripts_with_CDS.index))
    transcripts_with_CDS = pd.concat([transcripts_with_CDS, transcripts_cds_determined], ignore_index=True)
    transcripts_with_CDS['end ORF'] = transcripts_with_CDS['name']\
        .str.extract(r'[A-Z0-9\.]*:[A-Z0-9\.]*:ORF\-\d*:\d*\-(\d)')
    print("len together", len(transcripts_with_CDS.index))
    
    last_exon_length_dict = get_length_last_exon(transcripts_with_CDS['tid'].to_list(), gtf_file)
    print(len(last_exon_length_dict) == len(transcripts_with_CDS.index))
    transcripts_with_CDS['last_exon_length'] = transcripts_with_CDS['tid'].map(last_exon_length_dict)
    print(transcripts_with_CDS.head())


    pb.cleanup(remove_all=True)


    print('time needed for fasta', get_fasta)
    print('time needed for ORFs', getORFs)
    print('time needed for genomic', genomic_time)
    print('time needed for reference preparation', prepare_reference)
    print('time needed for coinc start sites', t_co_ss)
    print('Time needed to filter beds', t_filter_beds)
    print('Time needed to intersect and sum counts', time_intersect_sum)


    return transcripts_with_CDS

    