import sys

import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#import pybedtools as pb


from new_helper_functions import find_termination_codon, compose_transcript, get_transcript_string
from cds_determination_old import get_fasta_tid, get_cds_genomic_coordinates
from OrfFinder_py3_old import OrfFinder

def main():
    #read in gtf files
    #ensembl_gtf = GTF(sys.argv[2], check_ensembl_format=False)
    custom_gtf = GTF(sys.argv[1], check_ensembl_format=False)

    #get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    print("number of transcripts to investigate", len(transcript_ids))


    #build dataframe to store the computed features
    NMD_features_df = pd.DataFrame(columns = ["50_nt", "has_cds"],
                                    index = transcript_ids)

    #select information on transcripts as a list
    #transcript_string = get_transcript_string(transcript_ids)
    transcript_gtftk_object = custom_gtf\
        .select_by_key('feature', "CDS,exon,stop_codon")\
        .extract_data('transcript_id,start,end,exon_number,feature,strand,chrom,gene_id,score',
                       as_dict_of_merged_list=True)
    del custom_gtf
    transcript_ids_wo_cds = []

    for transcript_id in transcript_ids:
        transcript_info = transcript_gtftk_object[transcript_id]
        #build sublist for each feature of the transcript
        transcript_info = [transcript_info[x:x+8] for x in range(0, len(transcript_info), 8)]
        #extract cds features
        cds = [sub_list for sub_list in transcript_info if "CDS" in sub_list]
        #if there are CDS features: CDS is defined
        if len(cds) > 0:
            NMD_features_df.loc[transcript_id,"has_cds"] = 1
            stop_pos_genome = find_termination_codon(transcript_info, cds)
            stop_pos_transcript, last_ejc = compose_transcript(transcript_info, stop_pos_genome)
            if (last_ejc - stop_pos_transcript) >= 50:
                NMD_features_df.loc[transcript_id, "50_nt"] = 1
            else:
                NMD_features_df.loc[transcript_id, "50_nt"] = 0
        else: 
            NMD_features_df.loc[transcript_id,"has_cds"] = 0
            transcript_ids_wo_cds.append(transcript_id)
    
    #get fasta of transcripts with known id
    if len(transcript_ids_wo_cds) > 0:
        print("known_tids_no_cds", len(transcript_ids_wo_cds))
        transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys()\
                               if k in transcript_ids_wo_cds[0:1000]}
        genome_file = "./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa"
        sequences = get_fasta_tid(transcripts_no_cds, genome_file)
        ORFs = OrfFinder(sequences)
        cds_bed_positions, start_positions = get_cds_genomic_coordinates(ORFs)
        
        
        #extract gene ids of the transcripts that make ORFs
        gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
        gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
        #prepare string of gene ids to filter reference anno for gnees
        gene_string = get_transcript_string(gene_ids_ORF_transcripts)
        #load reference annotation
        ensembl_gtf = GTF(sys.argv[2], check_ensembl_format=False)
        #filter ensembl for gene_ids and protein coding transcripts
        reference_genes = ensembl_gtf\
        .select_by_key("gene_id", gene_string)\
        .select_by_key('feature', "CDS,exon,stop_codon,start_codon")\
        .extract_data('gene_id,transcript_id,start,end,exon_number,feature,strand,chrom,score',
                       as_dict_of_merged_list=True)
        del ensembl_gtf

        for start in start_positions.split("\n"):
            print(start)
            print(re.split(r"\t|:", start))


        for gene in gene_ids_ORF_transcripts[0:10]:
            try:
                gene_info = reference_genes[gene]
                #build sublist for each feature of the transcript
                gene_info = [gene_info[x:x+8] for x in range(0, len(gene_info), 8)]
                tids = [gene[0] for gene in gene_info]
                tids = list(set(tids))


                #separate information per transcript
                gene_dict = {}
                start_sites = []
                for tid in tids:
                    gene_dict[tid] = [gene for gene in gene_info if gene[0] == tid]
                    start_sites += [sub_list for sub_list in gene_dict[tid] if "start_codon" in sub_list]
                print(start_sites)
                #extract cds features
                #cds = [sub_list for sub_list in transcript_info if "CDS" in sub_list]
                #stop
                #start


            except KeyError as e:
                print(e, "gene", gene, "is not in the reference")


        #filter for protein coding transcripts with CDS, stop and start codon
        #check whether there are several of the same gene for which the start codon coincides
            #with the start of the ORF from the possible NMD transcript

        #bedtools intersect: select the ORF with the highest coverage in case that there are several left
        #how to deal with Ensmebl CDS annotation? Only take entries where they have a start and a stop annotated?
    #print(len(NMD_features_df["50_nt"]), sum(NMD_features_df["50_nt"]))
    return  0





if __name__ == "__main__":
    main()
