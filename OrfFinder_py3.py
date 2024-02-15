# Orfinder script written by Marcel Schulz 2019 adapted to py3.XX by Christina Kalk 2024
# Orffinder finds all possible ORFs of a minimal length in a given transcript and outputs
# a list of Sequence Records from the Biopyhton package is returned
#the header naem is composed of the following:
# id(current header of the transcript):ORF-number_of_ORF_in_the_current_transcript:start_position_in_transcript:stop_position_in_transcript

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



minOrfLength=50  #minimum number of amino acids per ORF

#original codons functions by natasha.sernova obtained from Biostars:
#https://www.biostars.org/p/229060/
#code has been modified

def get_orfs(seq, id, gene, description, annotations, countOrfs, orf_records):

        stops = ["TAA","TGA","TAG"]    
        lst1 = [] #List for the start codons
        lst2 = [] #List for the stop codons
        start = 0 #The start position of the sequence.
        counter = 0 #Counter for 3 optional orfs.
        #initializes the lists for 3 optional orfs.
        for i in range(3):
            lst1.append([])
            lst2.append([])
        #Add to the lists the positions of the start and stop codons.
        while (seq and counter < 3):

            for i in range(start, len(seq)-2, 3):
                codon = seq[i:i+3] #The codon is 3 nucleotides.
                #print codon+ "\t"
                if(codon == "ATG"): #If the codon is  a start codon.
                    lst1[start].append(i) #The position of the start codon.

                if(codon in stops): #if the codon is a stop codon.
                    lst2[start].append(i) #The position of the stop codon.


            start += 1 #promotes the starting position.
            counter += 1 #promotes the counter

        #for each reading frame go through the start site and extract and output possible proteins
        for frame in range(3):
                if len(lst1[frame]) > 0 and len(lst2[frame]) > 0:  #at least one start and stop codon per frame must exist
                    currentStart=0 #the current start position of a start codon in the frame
                    currentStop=0  #the current start position of a stop codon in the frame
                    while  currentStart < len(lst1[frame]) and currentStop < len(lst2[frame]): #codons are available
                        if lst1[frame][currentStart] < lst2[frame][currentStop] :  #found valid ORF
                                orf = seq[lst1[frame][currentStart]:lst2[frame][currentStop]] #this excludes the first stop position, but includes the start
                                if (len(orf)%3 == 0) and not ("N" in orf) and (len(orf)/3 >= minOrfLength):
                                    ##### added Frame to ID #####
                                    header = ''.join([id, ":ORF-", str(countOrfs), ":", str(lst1[frame][currentStart]), ":", str(lst2[frame][currentStop]+2)])
                                    sequence = seq[lst1[frame][currentStart]:lst2[frame][currentStop]+3]
                                    countOrfs = countOrfs + 1
                                    orf_records.append(SeqRecord(id= header, seq = sequence, name = gene, description = description, annotations = annotations))
                                #remove all other start codons that are nested between the current start and stop
                                while currentStart < len(lst1[frame]) and lst1[frame][currentStart] < lst2[frame][currentStop]:
                                        currentStart=currentStart+1
                                currentStop=currentStop+1
                                
                        elif lst1[frame][currentStart] > lst2[frame][currentStop]:
                                currentStop=currentStop+1 
        return(orf_records)        



def OrfFinder(transcript_SeqRecords):
    orf_records = []
    for transcript in transcript_SeqRecords:
            countOrfs=1#ORF numbers per transcript
            orf_records = get_orfs(transcript.seq.upper(), transcript.id, transcript.name, transcript.description, transcript.annotations, countOrfs, orf_records)
    return(orf_records)
                
        

