from pygtftk.gtf_interface import GTF
import sys



#select all transcript and exon entries from ensembl gtf and write to a new gtf file
def main():
    ensembl_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    transcript_ids = ensembl_gtf.select_by_key('feature', 'CDS').extract_data("transcript_id", as_list=True)
    transcript_string = ''
    for transcript in transcript_ids[:-1]:
        transcript_string += transcript+","
    transcript_string += transcript_ids[-1]
    no_nmd = ensembl_gtf.select_by_key('feature', 'transcript,exon').select_by_key('transcript_id', transcript_string)#.select_by_key('transcript_biotype', 'nonsense_mediated_decay', invert_match= True)
    
    no_nmd.write(sys.argv[2])
   

if __name__ == "__main__":
    main()




   