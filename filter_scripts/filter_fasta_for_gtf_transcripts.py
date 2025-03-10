import sys
from Bio import SeqIO
from pygtftk.gtf_interface import GTF

fasta_file = sys.argv[1]
gtf_file = sys.argv[2]
outfasta = sys.argv[3]
version_numbers = sys.argv[4]


def filter_fasta_by_gtf(gtf_file, fasta_file, outfasta, version_numbers=False):

    gtf = GTF(gtf_file, check_ensembl_format=False)
    gtf_tids = gtf.get_tx_ids(nr=True)

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    if version_numbers:
        records_to_keep = {trans_id: seq for (
            trans_id, seq) in record_dict.items() if trans_id.split('.')[0] in gtf_tids}

    else:
        records_to_keep = {trans_id: seq for (
            trans_id, seq) in record_dict.items() if trans_id in gtf_tids}

    with open(outfasta, "w") as output_handle:
        SeqIO.write(records_to_keep, output_handle, "fasta")
