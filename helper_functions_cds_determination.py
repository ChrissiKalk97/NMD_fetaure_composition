import regex as re
from operator import itemgetter
from typing import Dict, List
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pybedtools
from pybedtools import BedTool
from pygtftk.gtf_interface import GTF

def get_transcript_string(transcript_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    transcript_string = ''
    if len(transcript_ids) > 1: 
        for transcript_id in transcript_ids[:-1]:
            transcript_string += transcript_id+','
    transcript_string += transcript_ids[-1]
    return transcript_string



def get_fasta_tid(transcripts_no_cds, genome_file):
    '''creates the fasta sequence of a transcript from its exon coordinates 
    of the gtf file '''
    genome_dict = SeqIO.index(genome_file, 'fasta')
    sequences = []
    for transcript_id, transcript_info in transcripts_no_cds.items():
        fasta_string = ''
        description = ''
        transcript_info = [transcript_info[x:x+8] for x in range(0, len(transcript_info), 8)]
        transcript_exons = [sub_list for sub_list in transcript_info if 'exon' in sub_list]
        #sort list according to start positions of exons
        transcript_exons = sorted(transcript_exons, key = itemgetter(int(0)))
        for exon in transcript_exons:
            #build the transcript sequence from the exons in 5' to 3' order (+-strand), exon by exon
            fasta_string += genome_dict[exon[5]].seq[int(exon[0])-1:int(exon[1])]#get sequence of the exon by 
            #-1: 1-based system as in Ensembl, but string indexing is 0-based
            #might need to provide this to be changed, if assemblies have used different annotations, e.g. NCBI

            #subsetting the chromosome at the respecitve start and stop positions
            #add exon number and genomic start end to the description
            description += 'exon' + exon[2] + '-' + exon[0] + '-' + exon[1] + ':'
        #add chromosome and strand in front of description, remove the last ':'
        description = 'chrom' + exon[5] + ':' + 'strand' + exon[4] + ':' + description[:-1]
        if exon[4] == '-':
            fasta_string = Seq(fasta_string).reverse_complement()
        sequences.append(SeqRecord(id = transcript_id, seq = fasta_string, name = exon[6],\
                                    description = description, annotations={'score': exon[7]}))#name is gene
    return sequences
    

def get_cds_genomic_coordinates(orf_sequences):
    '''Calculates the genomic positions of the ORFs found in the ORF finder
    step and stores them in a bed object, the exons and partial exons 
    constituting the ORF are listed'''
    bed_string = ''
    start_positions = ''
    transcripts = []
    for orf in orf_sequences:
        transcripts.append(orf.id.split(':')[0])
        name = orf.name + ':' + orf.id
        orf_start = int(orf.id.split(':')[2])
        orf_end = int(orf.id.split(':')[3])
        transcript_info = orf.description.split(':')
        chromosome = re.match('chrom(.+)', transcript_info[0]).group(1)
        strand = re.match('strand(.+)', transcript_info[1]).group(1)
        length = 0
        start_counter = 0
        if strand == '+':
            counter = 2
            while orf_end > 0:
                exon = transcript_info[counter]
                exon_start = int(exon.split('-')[1])
                exon_end = int(exon.split('-')[2])
                former_length = length
                length += (exon_end - exon_start)+1#exon end and start are included
                if orf_start < length: 
                #this strictly smaller is needed here, otherwise the error occurs that the stop is smaller than the start
                    if start_counter == 0:
                            #note down the start position
                            start_positions += chromosome + '\t' + str(exon_start + orf_start - former_length) + '\t' +\
                            str(exon_start + orf_start - former_length + 2) + '\t' + name +'\t' + orf.annotations['score']\
                                  + '\t' + strand+ '\n'
                            start_counter += 1
                    if exon_start + orf_end - former_length > exon_end:
                        bed_string += chromosome + '\t' + str(exon_start + orf_start - former_length) + '\t' +\
                            str(exon_end) + '\t' + name +'\t' + orf.annotations['score'] + '\t' + strand + '\n'
                        orf_start = length#promote orf start such that it is former length 
                        #for the next exon and we start at the first exon position
                    else:
                        bed_string += chromosome + '\t' + str(exon_start + orf_start - former_length)\
                              + '\t' + str(exon_start + orf_end - former_length) + '\t' + name + '\t'\
                                  + orf.annotations['score'] + '\t' + strand + '\n'
                        orf_end = 0
                counter += 1
        else:
            counter = len(transcript_info) - 1
            while orf_end > 0:
                exon =  transcript_info[counter]
                exon_start = int(exon.split('-')[1])
                exon_end = int(exon.split('-')[2])
                former_length = length
                length += (exon_end - exon_start)+1
                if orf_start < length:
                    if start_counter == 0:
                        start_positions += chromosome + '\t' + str(exon_end - orf_start + former_length - 2)\
                              + '\t' + str(exon_end - orf_start + former_length) +\
                            '\t' + name +'\t' + orf.annotations['score'] + '\t' + strand + '\n'
                        start_counter += 1
                    if exon_end - orf_end + former_length < exon_start:
                        bed_string += chromosome + '\t' + str(exon_start) + '\t' +\
                              str(exon_end - orf_start + former_length) +\
                            '\t' + name +'\t' + orf.annotations['score'] + '\t' + strand + '\n'
                        orf_start = length
                    else:
                        bed_string += chromosome + '\t' + str(exon_end - orf_end + former_length) +\
                            '\t' + str(exon_end - orf_start + former_length) + '\t' + name + '\t' +\
                                  orf.annotations['score'] + '\t' + strand + '\n'
                        orf_end = 0
                counter -= 1
                
    try:
        bed_string = bed_string[:-1]
        start_positions = start_positions[:-1]
        tids_start = []
        tids_bed = []
        for start in start_positions.split('\n'):
            tids_start.append(re.split(r'\t|:', start)[4])
        
        for entry in bed_string.split('\n'):
            tids_bed.append(re.split(r'\t|:', entry)[4])
        
        print('nr of transcripts that were in ORF file', len(set(transcripts)))
        print('nr tids which are in bed', len(set(tids_bed)))
        print('nr tids for which start noted', len(set(tids_start)))

        genomic_coordinates = pybedtools.BedTool(bed_string, from_string=True)\
            .saveas('genomic_ccordinates_ORFs.bed')
        return genomic_coordinates, start_positions
    except Exception as e:
        print(e)
        for line in bed_string.split('\n'):
            try:
                if line.split('\t')[1] > line.split('\t')[2]:
                    print(line)
            except:
                print(line)




def get_pct_reference(ref_name, gene_ids_ORF_transcripts):
        #prepare string of gene ids to filter reference anno for gnees
        gene_string = get_transcript_string(gene_ids_ORF_transcripts)
        
        #load reference annotation
        #filter ensembl for gene_ids and protein coding transcripts
        reference_gtf = GTF(ref_name, check_ensembl_format=False)
        reference_genes = reference_gtf\
        .select_by_key('gene_id', gene_string)\
        .select_by_key('transcript_biotype', 'protein_coding')\
        .select_by_key('feature', 'CDS,start_codon')\
        .extract_data('gene_id,transcript_id,start,end,exon_number,feature,strand,chrom,score',
                       as_dict_of_merged_list=True)
        return reference_genes, reference_gtf

def get_ORF_start_by_gene(start_positions):
    start_sites_by_gene = {}
    for start in start_positions.split('\n'):
        start_info = re.split(r'\t|:', start)
        try:
            if start_info[3] not in start_sites_by_gene.keys():
                start_sites_by_gene[start_info[3]] = [start_info]
            else:
                start_sites_by_gene[start_info[3]].append(start_info)
        except:
            print(start)
    return start_sites_by_gene
    
def get_reference_gene_info(reference_genes: Dict, gene: str) -> tuple[List[List[str]], List[str]]:
    '''obtain all entries of reference gtf belonging to the gene considered and extract all transcript
    ids of the gene'''
    gene_info = reference_genes[gene]
    #build sublist for each feature of the transcript
    gene_info = [gene_info[x:x+8] for x in range(0, len(gene_info), 8)]
    tids = [gene[0] for gene in gene_info]
    tids = list(set(tids))
    return gene_info, tids

def get_cds_start(cds : List[List[str]], strand : str):
    '''obtain CDS start position for a transcript's CDS 
    information given as a list a lists'''
    start_position_plus = float('inf')
    start_position_minus = 0
    for partial_cds in cds:
        if strand == '+':
            five_prime = int(partial_cds[1])
            if start_position_plus < five_prime:
                start_position_plus = five_prime
                partial_cds[2] = str(int(partial_cds[1])+2)
        else:
            five_prime = int(partial_cds[2])
            if start_position_minus > five_prime:
                start_position_minus = five_prime
                partial_cds[1] = str(int(partial_cds[2])-2)
    return partial_cds


def get_reference_start_sites(gene_info, tids: List[str]) -> List[List[str]]:
    '''Calcualte the start sites of the CDS of the reference transcripts belonging
      to the gene considerd'''
    #separate information per transcript
    gene_dict = {}
    start_sites = []
    for tid in tids:
        gene_dict[tid] = [gene for gene in gene_info if gene[0] == tid]
        start_sites += [sub_list for sub_list in gene_dict[tid] if 'start_codon' in sub_list]
        if not start_sites:
            cds = [sub_list for sub_list in gene_dict[tid] if 'CDS' in sub_list]
            try:
                start_sites.append(get_cds_start(cds, cds[0][5]))
            except:
                print('No start sites and no CDS found!')
                #seems like this never happens when focusing on protein_coding transcripts
    return start_sites

def find_orfs_to_consider(start_sites: List[List[str]], orfs_starts_for_gene):
    '''Take start sites of reference and check which ORF start sites coincide: consider these ORFs
    Also consider ORFS for which no coinc start site if no other ORF with css found for the transcript'''
    starts = [start[1] for start in start_sites]
    orf_with_coinciding_start = [orf for orf in orfs_starts_for_gene if orf[1] in starts] 
    #print('start_sites', start_sites)
    #for orf in orf_with_coinciding_start:
        #print(orf)
        #print("matching starts", [start for start in start_sites if start[1] == orf[1]])
    orf_wo_coinciding_start = [orf for orf in orfs_starts_for_gene if orf[1] not in starts]
    tids_with_starts = set([start[4] for start in orf_with_coinciding_start])
    #add orfs with tids that are not present in the orfs with coinc start sites 
    orfs_to_consider =  orf_with_coinciding_start +\
            [orf for orf in orf_wo_coinciding_start if orf[4] not in tids_with_starts]
    #print("orfs to consider", orfs_to_consider)
    return orfs_to_consider

def get_transcript_orf_dict(orfs_to_consider):
    transcript_dict = {}
    for orf in orfs_to_consider:
        transcript_id = orf[4]
        if transcript_id not in transcript_dict.keys():
            transcript_dict[transcript_id] = [orf]
        else: 
            transcript_dict[transcript_id].append(orf)
    return transcript_dict

def categorize_transcripts_by_nr_orfs(transcript_dict: dict, transcripts_cds_determined: list,\
        transcripts_several_orfs: list, genes_tids_several_orfs: List[str], gene: str)\
            -> tuple[List[str], List[str], List[dict]]:
    for transcript, orfs in transcript_dict.items():
        if len(orfs) > 1:
            orf_list = [orf[3]+':'+orf[4]+':'+orf[5]+':'+orf[6]+':'+orf[7] for orf in orfs]
            transcripts_several_orfs = transcripts_several_orfs + orf_list
            genes_tids_several_orfs.append(gene)
        else:
            #write the single ORF as a dict
            transcripts_cds_determined.append({'tid': transcript,\
            'name': orfs[0][3] + ':' + orfs[0][4]+':' + orfs[0][5]\
                    + ':' + orfs[0][6] + ':' + orfs[0][7] , 'name_tar': ''})
    return transcripts_several_orfs, genes_tids_several_orfs, transcripts_cds_determined

def coinciding_start_sites(gene_ids_ORF_transcripts, reference_genes, orf_start_sites_by_gene):
    '''Determines CDS for transcripts where only one ORF is considered, notes down ORFs where several
    need to be considered, filters for transcripts without protein coding reference'''
    transcripts_cds_determined = []
    transcripts_several_orfs = []
    genes_tids_several_orfs = []
    for gene in gene_ids_ORF_transcripts:
        try:
            gene_info, tids = get_reference_gene_info(reference_genes, gene)
            
            #get predicted ORf start sites for the corresponding gene
            orfs_starts_for_gene = orf_start_sites_by_gene[gene]
            
            #get start sites of the protein coding reference CDS
            start_sites = get_reference_start_sites(gene_info, tids)
            if start_sites:
                orfs_to_consider = find_orfs_to_consider(start_sites, orfs_starts_for_gene)
                transcript_dict = get_transcript_orf_dict(orfs_to_consider)
                transcripts_several_orfs, genes_tids_several_orfs, transcripts_cds_determined\
                  =  categorize_transcripts_by_nr_orfs(transcript_dict, transcripts_cds_determined,\
                    transcripts_several_orfs, genes_tids_several_orfs, gene)
            
                #calculation of tids without coinciding start sites: not needed just info
                #tids_wo_starts = list(set([start[4] for start in orf_wo_coinciding_start]))
                #tids_orf_no_coinciding_start.extend([tid for tid in tids_wo_starts\
                                                      #if tid not in tids_with_starts])
        except KeyError as e:
            pass
    
    transcripts_cds_determined = pd.DataFrame(transcripts_cds_determined)
    return transcripts_cds_determined, transcripts_several_orfs, genes_tids_several_orfs



def filter_bed_file(transcriptIds, bedfile):
    '''Filter a bed file according to the name attribute being present in a set'''
    bed_of_transcripts = bedfile.filter(lambda gene:\
                                gene.name in set(transcriptIds))
    return bed_of_transcripts




    
def intersect_and_select_greatest_overlap(transcripts_several_orfs_bed, pc_reference_bed):#s: enforce strandedness
    intersection = transcripts_several_orfs_bed.intersect(pc_reference_bed, wao = True, s = True).saveas('intersection.bed')
    summed_counts = pd.read_table(intersection.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand',\
                    'chrom_tar', 'start_tar', 'stop_tar', 'name_tar', 'score_tar', 'strand_tar', 'overlap'], low_memory=False)
    
    summed_counts = summed_counts.groupby(['name', 'name_tar'])['overlap'].sum()
    summed_counts = summed_counts.reset_index()
    
    
    summed_counts['tid'] = summed_counts['name'].str.extract(r'[A-Z0-9\.]*:([A-Z0-9\.]*):.*')
    rowIds = summed_counts.groupby('tid')['overlap'].idxmax()
    transcripts_with_CDS = summed_counts.loc[rowIds]
    transcripts_with_CDS = transcripts_with_CDS.reset_index()[['name', 'name_tar', 'tid']]
    #append the others with only one CDS 
    return transcripts_with_CDS



def get_length_last_exon(transcript_ids_list, gtf_file):
    '''extract last exon length per transcript from reference'''
    gtf_file = GTF(gtf_file, check_ensembl_format=False)
    last_exon_length_by_transcript = {}
    transcript_string = get_transcript_string(transcript_ids_list)
    transcript_entries = gtf_file\
        .select_by_key('transcript_id', transcript_string)\
        .select_by_key('feature', 'exon')\
        .extract_data('transcript_id,start,end,exon_number,feature,strand',
                       as_dict_of_merged_list=True)
    for transcript_id in transcript_ids_list:
        transcript_exons = transcript_entries[transcript_id]
        transcript_exons = [transcript_exons[x:x+5] for x in range(0, len(transcript_exons), 5)]
        strand = transcript_exons[0][4]
        #sort transcript exons in ascending order
        transcript_exons = sorted(transcript_exons, key = itemgetter(int(0)))
        if strand == '+':
            last_exon = transcript_exons[-1]
        else:
            last_exon = transcript_exons[0]
        last_exon_length_by_transcript[transcript_id] = int(last_exon[1]) - int(last_exon[0])
    return last_exon_length_by_transcript

def calculate_50nt_rule(transcripts_with_CDS: pd.DataFrame, sequences) -> pd.DataFrame:
    t_length_dict = {}
    for sequence_rec in sequences:
        t_length_dict[sequence_rec.id] = len(sequence_rec.seq)
    
    transcripts_with_CDS['t_length'] = transcripts_with_CDS['tid'].map(t_length_dict)
    transcripts_with_CDS.set_index('tid', inplace = True)
    transcripts_with_CDS['distance_stop_EJC'] = \
    transcripts_with_CDS['t_length'].astype(int) -\
    transcripts_with_CDS['last_exon_length'].astype(int) -\
    transcripts_with_CDS['end_ORF'].astype(int)
    transcripts_with_CDS['50_nt_rule'] = transcripts_with_CDS['distance_stop_EJC'].ge(50).astype(int)
    return transcripts_with_CDS



                    
#there is the -s option for bedtools to enforce strandedness: overlaps are only reported if on the same strand