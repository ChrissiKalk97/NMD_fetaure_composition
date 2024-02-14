import regex as re
from operator import itemgetter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pybedtools
from pybedtools import BedTool
from pygtftk.gtf_interface import GTF

from new_helper_functions import get_transcript_string

def get_fasta_tid(transcripts_no_cds, genome_file):
    """creates the fasta sequence of a transcript from its exon coordinates 
    of the gtf file """
    genome_dict = SeqIO.index(genome_file, "fasta")
    sequences = []
    for transcript_id, transcript_info in transcripts_no_cds.items():
        fasta_string = ""
        description = ""
        transcript_info = [transcript_info[x:x+8] for x in range(0, len(transcript_info), 8)]
        transcript_exons = [sub_list for sub_list in transcript_info if "exon" in sub_list]
        transcript_exons = sorted(transcript_exons, key = itemgetter(int(0)))#sort list according to start positions of exons
        for exon in transcript_exons:
            #build the transcript sequence from the exons in 5' to 3' order (+-strand), exon by exon
            fasta_string += genome_dict[exon[5]].seq[int(exon[0])-1:int(exon[1])]#get sequence of the exon by 
            #-1: 1-based system as in Ensembl, but string indexing is 0-based
            #might need to provide this to be changed, if assemblies have used different annotations, e.g. NCBI

            #subsetting the chromosome at the respecitve start and stop positions
            description += "exon" + exon[2] + "-" + exon[0] + "-" + exon[1] + ":"#add exon number and genomic start end to the description
        description = "chrom" + exon[5] + ":" + "strand" + exon[4] + ":" + description[:-1]#add chromosome and strand in front of description,
        #remove the last ":"
        if exon[4] == "-":
            fasta_string = Seq(fasta_string).reverse_complement()
        sequences.append(SeqRecord(id = transcript_id, seq = fasta_string, name = exon[6], description = description, annotations={"score": exon[7]}))#name is gene
    return sequences
    

def get_cds_genomic_coordinates(orf_sequences):
    """Calculates the genomic positions of the ORFs found in the ORF finder step and stores them in a bed
    object, the exons and partial exons constituting the ORF are listed"""
    bed_string = ""
    start_positions = ""
    for orf in orf_sequences:
        name = orf.name + ":" + orf.id
        orf_start = int(orf.id.split(":")[2])
        orf_end = int(orf.id.split(":")[3])
        transcript_info = orf.description.split(":")
        chromosome = re.match("chrom(.+)", transcript_info[0]).group(1)
        strand = re.match("strand(.+)", transcript_info[1]).group(1)
        length = 0
        start_counter = 0
        if strand == "+":
            counter = 2
            while orf_end > 0:
                exon = transcript_info[counter]
                exon_start = int(exon.split("-")[1])
                exon_end = int(exon.split("-")[2])
                former_length = length
                length += (exon_end - exon_start)+1#exon end and start are included
                if orf_start < length: #this strictly smaller is needed here, otherwise the error occurs that the stop is smaller than the start
                    if start_counter == 0:
                            #note down the start position
                            start_positions += chromosome + "\t" + str(exon_start + orf_start - former_length) + "\t" +\
                            str(exon_start + orf_start - former_length + 2) + "\t" + name +"\t" + orf.annotations["score"] + "\t" + strand+ "\n"
                            start_counter += 1
                    if exon_start + orf_end - former_length > exon_end:
                        bed_string += chromosome + "\t" + str(exon_start + orf_start - former_length) + "\t" +\
                            str(exon_end) + "\t" + name +"\t" + orf.annotations["score"] + "\t" + strand + "\n"
                        orf_start = length#promote orf start such that it is former length for the next exon and we start at the first exon position
                    else:
                        bed_string += chromosome + "\t" + str(exon_start + orf_start - former_length) + "\t" +\
                            str(exon_start + orf_end - former_length) + "\t" + name +"\t" + orf.annotations["score"] + "\t" + strand + "\n"
                        orf_end = 0
                counter += 1
        else:
            counter = len(transcript_info) - 1
            while orf_end > 0:
                exon =  transcript_info[counter]
                exon_start = int(exon.split("-")[1])
                exon_end = int(exon.split("-")[2])
                former_length = length
                length += (exon_end - exon_start)+1
                if orf_start < length:
                    if start_counter == 0:
                        start_positions += chromosome + "\t" + str(exon_end - orf_start + former_length - 2) + "\t" + str(exon_end - orf_start + former_length) +\
                            "\t" + name +"\t" + orf.annotations["score"] + "\t" + strand + "\n"
                        start_counter += 1
                    if exon_end - orf_end + former_length < exon_start:
                        bed_string += chromosome + "\t" + str(exon_start) + "\t" + str(exon_end - orf_start + former_length) +\
                            "\t" + name +"\t" + orf.annotations["score"] + "\t" + strand + "\n"
                        orf_start = length
                    else:
                        bed_string += chromosome + "\t" + str(exon_end - orf_end + former_length) +\
                            "\t" + str(exon_end - orf_start + former_length) + "\t" + name + "\t" + orf.annotations["score"] + "\t" + strand + "\n"
                        orf_end = 0
                counter -= 1
                
    try:
        genomic_coordinates = pybedtools.BedTool(bed_string, from_string=True)
        #start_positions = pybedtools.BedTool(start_positions, from_string=True)
        return genomic_coordinates, start_positions
    except Exception as e:
        print(e)
        for line in bed_string.split("\n"):
            try:
                if line.split("\t")[1] > line.split("\t")[2]:
                    print(line)
            except:
                print(line)


def get_pct_reference(ref_name, gene_ids_ORF_transcripts):
        #prepare string of gene ids to filter reference anno for gnees
        gene_string = get_transcript_string(gene_ids_ORF_transcripts)
        #load reference annotation
        ensembl_gtf = GTF(ref_name, check_ensembl_format=False)
        #filter ensembl for gene_ids and protein coding transcripts
        reference_genes = ensembl_gtf\
        .select_by_key("gene_id", gene_string)\
        .select_by_key('feature', "CDS,exon,start_codon")\
        .extract_data('gene_id,transcript_id,start,end,exon_number,feature,strand,chrom,score',
                       as_dict_of_merged_list=True)
        del ensembl_gtf

        return reference_genes
    

                    
#there is the -s option for bedtools to enforce strandedness: overlaps are only reported if on the same strand