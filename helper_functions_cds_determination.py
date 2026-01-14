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
from cdifflib import CSequenceMatcher
import difflib  # import SequenceMatcher
difflib.SequenceMatcher = CSequenceMatcher


def get_transcript_string(transcript_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    transcript_string = ''
    if len(transcript_ids) > 1:
        for transcript_id in transcript_ids[:-1]:
            transcript_string += transcript_id+','
    transcript_string += transcript_ids[-1]
    return transcript_string


def get_fasta_tid(transcripts_no_cds, genome_file, seq_type: str, plus_stop=False):
    '''creates the fasta sequence of a transcript from its exon coordinates 
    of the gtf file '''
    genome_dict = SeqIO.index(genome_file, 'fasta')
    sequences = []
    for transcript_id, transcript_info in transcripts_no_cds.items():
        fasta_string = ''
        description = ''
        transcript_info = [transcript_info[x:x+8]
                           for x in range(0, len(transcript_info), 8)]
        transcript_exons = [
            sub_list for sub_list in transcript_info if seq_type in sub_list]
        # sort list according to start positions of exons
        # if len(transcript_exons) > 0:
        transcript_exons.sort(key=lambda elem: int(elem[0]))
        print(transcript_exons)
        for idx, exon in enumerate(transcript_exons):
            if plus_stop == True and idx == len(transcript_exons) - 1:
                fasta_string += genome_dict[exon[5]
                                            ].seq[int(exon[0])-1:int(exon[1])+3]
            else:
                # build the transcript sequence from the exons in 5' to 3' order (+-strand), exon by exon
                # get sequence of the exon by
                fasta_string += genome_dict[exon[5]
                                            ].seq[int(exon[0])-1:int(exon[1])]
                # -1: 1-based system as in Ensembl, but string indexing is 0-based
                # might need to provide this to be changed, if assemblies have used different annotations, e.g. NCBI

            # subsetting the chromosome at the respecitve start and stop positions
            # add exon number and genomic start end to the description
            description += 'exon' + exon[2] + \
                '-' + exon[0] + '-' + exon[1] + ':'
        # add chromosome and strand in front of description, remove the last ':'
        description = 'chrom' + exon[5] + ':' + \
            'strand' + exon[4] + ':' + description[:-1]
        if exon[4] == '-':
            fasta_string = Seq(fasta_string).reverse_complement()
        sequences.append(SeqRecord(id=transcript_id, seq=Seq(fasta_string), name=exon[6],
                                   # name is gene
                                   description=description, annotations={'score': exon[7]}))
    del genome_dict
    return sequences


def get_cds_genomic_coordinates(orf_sequences):
    '''Calculates the genomic positions of the ORFs found in the ORF finder
    step and stores them in a bed object, the exons and partial exons 
    constituting the ORF are listed'''
    bed_string = ''
    transcripts = []
    orf_dict_exon_with_stop_length = {}
    distance_stop_next_EJC = {}
    for orf in orf_sequences:
        tid = orf.id.split(':')[0]
        transcripts.append(tid)
        name = orf.name + '|' + orf.id
        orf_start = int(orf.id.split(':')[2])
        orf_end = int(orf.id.split(':')[3])
        transcript_info = orf.description.split(':')
        chromosome = re.match('chrom(.+)', transcript_info[0]).group(1)
        strand = re.match('strand(.+)', transcript_info[1]).group(1)
        length = 0
        if strand == '+':
            counter = 2
            while orf_end > 0:
                exon = transcript_info[counter]
                exon_start = int(exon.split('-')[1])
                exon_end = int(exon.split('-')[2])
                former_length = length
                # exon end and start are included
                length += (exon_end - exon_start)+1
                if orf_start < length:
                    # this strictly smaller is needed here, otherwise the error occurs that the stop is smaller than the start
                    if exon_start + orf_end - former_length > exon_end:
                        bed_string += chromosome + '\t' + str(exon_start + orf_start - former_length) + '\t' +\
                            str(exon_end) + '\t' + name + '\t' + \
                            orf.annotations['score'] + '\t' + strand + '\n'
                        orf_start = length  # promote orf start such that it is former length
                        # for the next exon and we start at the first exon position
                    else:
                        bed_string += chromosome + '\t' + str(exon_start + orf_start - former_length)\
                            + '\t' + str(exon_start + orf_end - former_length) + '\t' + name + '\t'\
                            + orf.annotations['score'] + '\t' + strand + '\n'
                        distance_stop_next_EJC[orf.name +
                                               '|' + orf.id] = length - orf_end + 1
                        orf_end = 0
                counter += 1
                # the last exon encountered has the stop codon, note down exon length
                orf_dict_exon_with_stop_length[orf.name +
                                               '|' + orf.id] = exon_end - exon_start + 1

        else:
            counter = len(transcript_info) - 1
            while orf_end > 0:
                exon = transcript_info[counter]
                exon_start = int(exon.split('-')[1])
                exon_end = int(exon.split('-')[2])
                former_length = length
                length += (exon_end - exon_start)+1
                if orf_start < length:
                    if exon_end - orf_end + former_length < exon_start:
                        bed_string += chromosome + '\t' + str(exon_start) + '\t' +\
                            str(exon_end - orf_start + former_length) +\
                            '\t' + name + '\t' + \
                            orf.annotations['score'] + '\t' + strand + '\n'
                        orf_start = length
                    else:
                        bed_string += chromosome + '\t' + str(exon_end - orf_end + former_length) +\
                            '\t' + str(exon_end - orf_start + former_length) + '\t' + name + '\t' +\
                            orf.annotations['score'] + '\t' + strand + '\n'
                        distance_stop_next_EJC[orf.name +
                                               '|' + orf.id] = length - orf_end + 1
                        orf_end = 0
                counter -= 1
                # the last exon encountered has the stop codon, note down exon length
                orf_dict_exon_with_stop_length[orf.name +
                                               '|' + orf.id] = exon_end - exon_start + 1

    try:
        bed_string = bed_string[:-1]
        tids_bed = []
        for entry in bed_string.split('\n'):
            tids_bed.append(re.split(r'\t|:|\|', entry)[4])

        print('nr of transcripts that were in ORF file', len(set(transcripts)))
        print('nr tids which are in bed', len(set(tids_bed)))

        genomic_coordinates = pybedtools.BedTool(bed_string, from_string=True)\
            .saveas('genomic_coordinates_ORFs.bed')
        return genomic_coordinates, orf_dict_exon_with_stop_length, distance_stop_next_EJC
    except Exception as e:
        print(e)
        for line in bed_string.split('\n'):
            try:
                if line.split('\t')[1] > line.split('\t')[2]:
                    print(line)
            except:
                print(line)


def get_pct_reference(ref_name, gene_ids_ORF_transcripts):
    # prepare string of gene ids to filter reference anno for gnees
    gene_string = get_transcript_string(gene_ids_ORF_transcripts)
    # load reference annotation
    # filter ensembl for gene_ids and protein coding transcripts
    # take CDS for the intersection as we require that the
    # CDs of the source transcripts overlaps with the reference
    # CDS
    reference_gtf = GTF(ref_name, check_ensembl_format=False)

    reference_gtf_CDS = reference_gtf\
        .select_by_key('gene_id', gene_string)\
        .select_by_key('transcript_biotype', 'protein_coding')\
        .select_by_key('feature', 'CDS')
    return reference_gtf_CDS


def longest_common_substring(source: str, target: str) -> str:
    """Computes the longest common substring of s1 and s2"""
    seq_matcher = difflib.SequenceMatcher(
        isjunk=None, a=source, b=target, autojunk=False)
    match = seq_matcher.find_longest_match(0, len(source), 0, len(target))

    if match.size:
        return source[match.a: match.a + match.size]
    else:
        return ""


def longest_common_substring_percentage(source: str, target: str) -> float:
    """Computes the longest common substring percentage of s1 and s2"""
    assert min(len(source), len(target)
               ) > 0, "One of the given string is empty"
    # get target exact match percentage
    longest_match = longest_common_substring(source, target)
    return (len(longest_match),
            # divide common substring by source length
            len(longest_match)/len(source))


def select_row(group):
    # Check if there's any entry above a certain value in the value_column
    high_overlap = group[group['protein_overlap_perc'] > 0.1]
    high_overlap = high_overlap[high_overlap['protein_overlap_aa'] >= 20]
    if not high_overlap.empty:
        high_overlap = high_overlap.sort_values(
            by=['protein_overlap_aa'], ascending=[False])
        # returns lowest start posiiton, but if tie
        return high_overlap.loc[high_overlap['start_ORF'].idxmin()]
    # sorting according to protein overlap decides
    else:
        # sort by AA overlap, such that if several ORFs start at the same position the longer one
        # is chosen??? Should not happen...
        group = group.sort_values(by=['protein_overlap_aa'], ascending=[False])
        return group.loc[group['start_ORF'].idxmin()]


def find_cds_orf(reference_gtf, orf_bed_positions, orf_file, transcript_file, use_fasta,  prot_cod_bed):
    ''''''
    if use_fasta == False:
        reference_bed = BedTool(reference_gtf.to_bed(
            name=('gene_id', 'transcript_id'), sep='|')).saveas('pc_reference.bed')
        intersection = orf_bed_positions.intersect(
            reference_bed, wao=True, s=True).saveas('intersection.bed')
    else:
        reference_bed = BedTool(prot_cod_bed)
        intersection = orf_bed_positions.intersect(
            reference_bed, wao=True).saveas('intersection.bed')
        # leave out the strand for now as the coordinates differ (+/- vs /-1)
    summed_overlap = pd.read_table(intersection.fn,
                                   names=['chrom', 'start', 'stop', 'name', 'score', 'strand',
                                          'chrom_tar', 'start_tar', 'stop_tar', 'name_tar', 'score_tar',
                                          'strand_tar', 'overlap'], low_memory=False)
    if use_fasta == True:
        summed_overlap['name_tar'] = summed_overlap['name_tar'] + \
            '|' + summed_overlap['score_tar']

    summed_overlap = summed_overlap.groupby(
        ['name', 'name_tar'])['overlap'].sum()
    summed_overlap = summed_overlap.reset_index()

    summed_overlap['tid'] = summed_overlap['name'].str.extract(
        r'[^|]*\|([^:]*):.*')

    # summed_overlap['tid'] = summed_overlap['name'].str.extract(
    #       r'[A-Z0-9\.]*\|([A-Z0-9\.]*):.*')

    print('before any filtering we have', len(
        summed_overlap['tid'].unique()), 'tids')

    # get gene ids of source and target, they need to be the same, otherwise the match is invalid
    # summed_overlap['gid'] = summed_overlap['name'].str.extract(
    #    r'([A-Z0-9\.]*)\|[A-Z0-9\.]*:.*')
    # summed_overlap['gid_target'] = summed_overlap['name_tar'].str.extract(
    #    r'([A-Z0-9\.]*)\|[A-Z0-9\.]*')
    summed_overlap['gid'] = summed_overlap['name'].str.extract(
        r'([^|]*)\|[^:]*:.*')
    summed_overlap['gid_target'] = summed_overlap['name_tar'].str.extract(
        r'([^|]*)\|.*')

    summed_overlap['ORF_nr'] = summed_overlap['name'].str.split(':').str[1]
    summed_overlap['start_ORF'] = summed_overlap['name'].str.split(':').str[2]
    summed_overlap['end_ORF'] = summed_overlap['name'].str.split(':').str[3]

    summed_overlap['start_ORF'] = summed_overlap['start_ORF'].astype(int)
    # end ORF is the first position of the stop codon
    summed_overlap['end_ORF'] = summed_overlap['end_ORF'].astype(int) - 2

    # get length of target transcripts
    """target_length_dict = {}
    for interval in reference_bed:
        if interval.name in target_length_dict.keys():
            target_length_dict[interval.name] += interval.stop - interval.start
        else:
            target_length_dict[interval.name] = interval.stop - interval.start
    summed_overlap['target_length'] = summed_overlap['name_tar'].map(
        target_length_dict)
    summed_overlap['target_coverage_percentage'] = summed_overlap['overlap'] / \
        summed_overlap['target_length']"""

    ORF_sequences = SeqIO.to_dict(SeqIO.parse(orf_file, "fasta"))
    target_sequences = SeqIO.to_dict(SeqIO.parse(transcript_file, "fasta"))

    summed_overlap['protein_overlap_perc'] = 0.0
    summed_overlap['protein_overlap_aa'] = 0
    for row in summed_overlap.itertuples():
        source = summed_overlap['name'].at[row.Index]
        target = summed_overlap['name_tar'].at[row.Index]
        if target != '.' and target != '.|.':
            source_seq = ORF_sequences[source]
            target_seq = target_sequences[target]
            longest_common_substring = \
                longest_common_substring_percentage(
                    str(source_seq.seq), str(target_seq.seq))
            summed_overlap['protein_overlap_perc'].at[row.Index]\
                = longest_common_substring[1]
            summed_overlap['protein_overlap_aa'].at[row.Index]\
                = longest_common_substring[0]

    # Apply the custom function to each group and concatenate the results
    transcripts_with_CDS = summed_overlap.groupby('tid').apply(select_row)

    return transcripts_with_CDS


# there is the -s option for bedtools to enforce strandedness: overlaps are only reported if on the same strand
