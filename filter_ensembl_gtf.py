from pygtftk.gtf_interface import GTF
import sys



#select all transcript and exon entries from ensembl gtf and write to a new gtf file
def main():
    ensembl_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    nmd_gtf = ensembl_gtf.select_by_key('feature', 'transcript,exon').select_by_key('transcript_biotype', 'nonsense_mediated_decay')
    nmd_gtf.write(sys.argv[2])
   

if __name__ == "__main__":
    main()




   