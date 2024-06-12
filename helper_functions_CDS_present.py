from typing import List
from operator import itemgetter

def check_stop_codon(transcript_as_list: List[List[str]], strand : str):
    stop_codon = [sub_list for sub_list in transcript_as_list if 'stop_codon' in sub_list]
    if len(stop_codon) > 0:
        if strand == '+':
            stop_codon = int(stop_codon[0][0])
        else:
            stop_codon = int(stop_codon[0][1])
        return stop_codon
    else:
        return None
    
def get_cds_end(cds : List[List[str]], strand : str)-> int:
    stop_position_plus = 0
    stop_position_minus = float('inf')
    cds_length = 0
    for partial_cds in cds:
        if strand == '+':
            three_prime = int(partial_cds[1])
            five_prime = int(partial_cds[0])
            if stop_position_plus < three_prime:
                stop_position_plus = three_prime
                cds_end = three_prime
                cds_length += three_prime - five_prime + 1
        else:
            three_prime = int(partial_cds[0])
            five_prime = int(partial_cds[1])
            if stop_position_minus > three_prime:
                stop_position_minus = three_prime
                cds_end = three_prime
                cds_length += five_prime - three_prime + 1
    return cds_end, cds_length
                
 
def find_termination_codon(transcript_as_list: List[str], cds: List[List[str]]) -> int:
    #get strand information and stop codon if annotated
    strand = transcript_as_list[0][4]
    
    #if stop codon annotated: return its first base position
    stop_position = check_stop_codon(transcript_as_list, strand)
    if stop_position is not None:
        cds_end, cds_length = get_cds_end(cds, strand)
        return stop_position, cds_length
    
    #if stop codon not annotated, compose CDS and extract end poisition of CDS
    else:
        cds_end, cds_length = get_cds_end(cds, strand)
        #get exons of transcript
        exons = [sub_list for sub_list in transcript_as_list if 'exon' in sub_list]
        exons = sorted(exons, key = itemgetter(int(0)))#sort list according to start positions of exons
        if strand == '+':
            #check whether the end of the CDS coincides with the end of an exon
            exon_end_is_cds_end = [exon for exon in exons if str(cds_end) == exon[1]]
            #if not assign baseposition on genome after cds (CDs_end+1) as stop
            if len(exon_end_is_cds_end) == 0:
                stop_position = cds_end+1#stop position first base will start right after the CDS
            #if CDS end coincides with exon end
            else:
                #extract exon number of coinciding exon
                exon_num = exons.index(exon_end_is_cds_end[0])
                #annotate the start of the next exon as stop position
                try:
                    next_exon = exons[exon_num+1]
                    #if cds ends with exon and there is a next exon, set stop position to first bp of next exon
                    stop_position = int(next_exon[0])
                except IndexError:
                    #stop_position = cds_end-2
                    stop_position = cds_end
        else:
            exon_end_is_cds_end = [exon for exon in exons if str(cds_end) == exon[0]]
            if len(exon_end_is_cds_end) == 0:
                stop_position = cds_end-1
            else:
                #extract exon number of coinciding exon
                exon_num = exons.index(exon_end_is_cds_end[0])
                #annotate the start of the next exon as stop position
                try:
                    if exon_num != 0:
                        next_exon = exons[exon_num-1]
                    else:
                        next_exon = exon_end_is_cds_end[0]
                    #if cds ends with exon and there is a next exon, set stop position to first bp of next exon
                    stop_position = int(next_exon[1])
                except IndexError:
                    #stop_position = cds_end+2
                    stop_position = cds_end
        return stop_position, cds_length
    

def compose_transcript(transcript_as_list, stop_position):
    #extract all exons of the transcript
    exons = [sub_list for sub_list in transcript_as_list if 'exon' in sub_list]
    #print(exons)
    strand = exons[0][4]
    #dict to store transcript coordinates for exon junctions and stop_position
    transcript_coordinates = {'length': 0, 'exon_jc': {}, 'stop_position': 0}
    last_exon = 0
    if strand == '+':
        exons.sort(key=lambda elem: int(elem[0]))#sort list according to start positions of exons
    else:
        exons.sort(key=lambda elem: int(elem[0]), reverse = True)
    exon_index = 0
    for exon in exons:
        exon_num = exon[2]
        if strand == '+':
            five_prime = int(exon[0])
            three_prime = int(exon[1])
            length = three_prime - five_prime +1
            #if the stop position is inside this exon: calculate the stop position in transcript coordinates, length plus distance from stop to 5' end of exon
            if stop_position >= five_prime and stop_position <= three_prime:
                transcript_coordinates['stop_position'] = transcript_coordinates['length'] + stop_position - five_prime
                exon_containing_stop = exon_index
                dist_stop_next_EJC = three_prime - stop_position

            transcript_coordinates['length']+= length#add length of exon to total length of the transcript
            transcript_coordinates['exon_jc'][exon_num] = length#note down the length of the transcript including current exon = EJC position
            if last_exon < int(exon_num):
                last_exon = int(exon_num)
    
        else:
            five_prime = int(exon[1])
            three_prime = int(exon[0])
            #for the negative strand the ordering is adjusted such that the first exon is the on with
            #the largest genoic +-strand cooridnates
            #the five prime end is bigger than the 3 prime end (as this is in respect to the negative strand)
            #plus 1 because the las coordinate is also part of the exon
            length = five_prime - three_prime + 1
            if stop_position <= five_prime and stop_position >= three_prime:
                #as cooridnates for stop and five_prime are on +strand but we are calculating on the -strand
                #distance from current (minus strand) five prime to PTC, is equivaluent to 5#-stopp_position in +strand coordinates
                transcript_coordinates['stop_position'] =  transcript_coordinates['length'] - stop_position + five_prime
                exon_containing_stop = exon_index
                dist_stop_next_EJC = five_prime - stop_position

            transcript_coordinates['length']+= length
            transcript_coordinates['exon_jc'][exon_num] = length
            if last_exon < int(exon_num):
                last_exon = int(exon_num)
        exon_index += 1
        
    
    trans_coord_last_ejc = transcript_coordinates['length'] - transcript_coordinates['exon_jc'][str(last_exon)]

    exon_containing_stop_length = int(exons[exon_containing_stop][1]) - int(exons[exon_containing_stop][0]) +1
    
    return transcript_coordinates, trans_coord_last_ejc, exon_containing_stop_length, dist_stop_next_EJC