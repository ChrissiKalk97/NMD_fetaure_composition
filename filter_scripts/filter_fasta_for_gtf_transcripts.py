import sys
from Bio import SeqIO
from pygtftk.gtf_interface import GTF


def filter_fasta_by_gtf(gtf_file, fasta_file, outfasta, version_numbers=False):

    gtf_lines = GTF(gtf_file, check_ensembl_format=False)
    gtf_tids = gtf_lines.get_tx_ids(nr=True)

    records_to_keep = []
    with open(fasta_file) as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if version_numbers and record.id.split('.')[0] in gtf_tids:
                records_to_keep.append(record)
            elif record.id in gtf_tids:
                records_to_keep.append(record)
    with open(outfasta, "w") as output_handle:
        SeqIO.write(records_to_keep, output_handle, "fasta")


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    gtf_file = sys.argv[2]
    outfasta = sys.argv[3]
    version_numbers = sys.argv[4]
    filter_fasta_by_gtf(gtf_file, fasta_file, outfasta, version_numbers)
