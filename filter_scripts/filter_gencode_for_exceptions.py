from pygtftk.gtf_interface import GTF
import sys



#select all transcript and exon entries from ensembl gtf and write to a new gtf file
def main():
    gencode_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    nmd_exceptions = gencode_gtf.select_by_key('feature', 'transcript').select_by_key('tag', 'NMD_exception')
    #.select_by_key(tag, 'mRNA_end_NF,cds_end_NF,mRNA_start_NF,cds_start_NF', invert_match = True)


    nmd_ex_transcripts = nmd_exceptions.extract_data("transcript_id", as_list=True)
    print('number of NMD exceptions', len(nmd_ex_transcripts))
    
    nmd_exceptions.write(sys.argv[2])
   

if __name__ == "__main__":
    main()




   