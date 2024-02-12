import sys


import pandas as pd
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#import pybedtools as pb

from new_helper_functions import find_termination_codon, compose_transcript
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
        transcripts_no_cds = {k: transcript_gtftk_object[k] for k in transcript_gtftk_object.keys() if k in transcript_ids_wo_cds[0:10000]}
        genome_file = "./Homo_sapiens.GRCh38.dna.primary_assembly_110_new.fa"
        sequences = get_fasta_tid(transcripts_no_cds, genome_file)
        ORFs = OrfFinder(sequences)
        cds_bed_positions = get_cds_genomic_coordinates(ORFs)
        #extract gene ids of the transcripts that make ORFs
        #gene_ids_ORF_transcripts = [orf.name for orf in ORFs]
        #gene_ids_ORF_transcripts = list(set(gene_ids_ORF_transcripts))
        #print(gene_ids_ORF_transcripts[0:100])
        #print(len(gene_ids_ORF_transcripts))
        
        #load ensembl gtf
        #ensembl_gtf = GTF(sys.argv[2], check_ensembl_format=False)
        #filter ensembl for gene_ids and protein coding transcripts: extract cds: plus and minus aware
        #check for coinciding start sites: filter out the ORFs for which no start site coincides
        #bedtools intersect: select the ORF with the highest coverage in case that there are several left
        #how to deal with Ensmebl CDS annotation? Only take entries where they have a start and a stop annotated?
    #print(len(NMD_features_df["50_nt"]), sum(NMD_features_df["50_nt"]))
    return  0





if __name__ == "__main__":
    main()
