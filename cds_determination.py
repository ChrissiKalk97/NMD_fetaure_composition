from helper_functions_cds_determination import get_fasta_tid, get_cds_genomic_coordinates,\
      get_pct_reference, get_ORF_start_by_gene
from new_helper_functions import get_cds_start
from OrfFinder_py3 import OrfFinder

def determine_cds(transcript_gtftk_object, transcript_ids_wo_cds, reference_file, reference_sequences_file):
    transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds[0:10000]}
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
    #print(start_sites_by_gene)

    no_starts_counter = 0     
    for gene in gene_ids_ORF_transcripts:
        try:
            gene_info = reference_genes[gene]
            orfs_starts_for_gene = orf_start_sites_by_gene[gene]
            #build sublist for each feature of the transcript
            gene_info = [gene_info[x:x+8] for x in range(0, len(gene_info), 8)]
            tids = [gene[0] for gene in gene_info]
            tids = list(set(tids))


            #separate information per transcript
            gene_dict = {}
            start_sites = []
            for tid in tids:
                gene_dict[tid] = [gene for gene in gene_info if gene[0] == tid]
                #print(gene_dict[tid])
                start_sites += [sub_list for sub_list in gene_dict[tid] if "start_codon" in sub_list]
                if not start_sites:
                    cds += [sub_list for sub_list in gene_dict[tid] if "CDS" in sub_list]
                    if cds:
                        start_sites.append(get_cds_start(cds, cds[0][5]))

            #print(orfs_starts_for_gene) 
            if start_sites:
                starts = [start[1] for start in start_sites]
                orf_with_coinciding_start = [orf for orf in orfs_starts_for_gene if orf[1] in starts] 
                orf_wo_coinciding_start = [orf for orf in orfs_starts_for_gene if orf[1] not in starts]   
                #print("with starts", len(orf_with_coinciding_start))
                #print("no starts", len(orf_wo_coinciding_start))
                #print("the starts", len(starts))
                tids_with_starts = list(set([start[4] for start in orf_with_coinciding_start]))
                tids_wo_starts = list(set([start[4] for start in orf_wo_coinciding_start]))
                no_starts_counter += len([tid for tid in tids_wo_starts if tid not in tids_with_starts])
            #no need to check the chromosome as we are on the same gene
            #print(f"transcript annotated CDS/start sites for {gene}", start_sites)
            #print(f"ORF starts for {gene}", orfs_starts_for_gene)

            
            
        except KeyError as e:
            print(e, "gene", gene, "is not in the reference")
        
        print("nr of transcripts with no starts", no_starts_counter)