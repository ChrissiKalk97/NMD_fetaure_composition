import sys
import pandas as pd
from pygtftk.gtf_interface import GTF

from new_helper_functions import get_transcript_string

#from Bio import SeqIO
#from cds_determination import get_transcript_fasta

#implementation for one transcript
def main():
    #read in gtf files
    ensembl_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    custom_gtf = GTF(sys.argv[2], check_ensembl_format=False)

    #get transcript ids
    transcript_ids = custom_gtf.get_tx_ids(nr=True)
    print("number of transcripts to investigate", len(transcript_ids))


    #select information on transcripts as a list
    transcript_string = get_transcript_string(transcript_ids)
    transcript_gtftk_object = ensembl_gtf.select_by_key('transcript_id', transcript_string)\
        .select_by_key('feature', "CDS,stop_codon,transcript")\
        .extract_data('transcript_id,start,end,exon_number,feature,strand',
                       as_dict_of_merged_list=True)

    CDS_ends_at_transcript = 0

    for transcript_id in transcript_ids:
        transcript_info = transcript_gtftk_object[transcript_id]
        #build sublist for each feature of the transcript
        transcript_info = [transcript_info[x:x+5] for x in range(0, len(transcript_info), 5)]
        #extract cds features
        stop_codon = [sub_list for sub_list in transcript_info if "stop_codon" in sub_list]
        cds = [sub_list for sub_list in transcript_info if "CDS" in sub_list]
        transcript = [sub_list for sub_list in transcript_info if "transcript" in sub_list]
        #print(transcript)
        if len(stop_codon) == 0:
            strand = transcript[0][4]
            for partial_cds in cds:
                stop_position_plus = 0
                stop_position_minus = float("inf")
                for partial_cds in cds:
                    if strand == "+":
                        three_prime = int(partial_cds[1])
                        if stop_position_plus < three_prime:
                            cds_end = three_prime
                    else:
                        three_prime = int(partial_cds[0])
                        if stop_position_minus > three_prime:
                            cds_end = three_prime
            if strand == "+":
                
                if cds_end == int(transcript[0][1]):
                    
                    CDS_ends_at_transcript += 1
            else:
                if cds_end == int(transcript[0][0]):
                    print(transcript_id, cds_end, transcript[0][1])
                    CDS_ends_at_transcript += 1
       
    print(CDS_ends_at_transcript, " transcripts were found where CDS end coincides with transcript end")
    return  0





if __name__ == "__main__":
    main()
