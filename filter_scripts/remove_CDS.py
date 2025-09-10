from pygtftk.gtf_interface import GTF
import sys


# select all transcript and exon entries from ensembl gtf and write to a new gtf file
def main():
    gtf = GTF(sys.argv[1], check_ensembl_format=False)
    filtered_gtf = gtf.select_by_key(
        'feature', 'CDS', invert_match=True)
    # .select_by_key(tag, 'mRNA_end_NF,cds_end_NF,mRNA_start_NF,cds_start_NF', invert_match = True)

    filtered_gtf.write(sys.argv[2])


if __name__ == "__main__":
    main()
