#This script takes a gtf of an assembly from stringtie that was merged of single sample assemblies
#and writes the ref_gene_id attribute as the gene_id
import sys
import regex as re

assembled_gtf = open(sys.argv[1], 'r')
Lines = assembled_gtf.readlines()
gtf_with_repalcement = open(sys.argv[2], 'w')

replaced_lines = []
count = 0

for line in Lines:
    ref_gene_id = re.search(r'ref_gene_id "(ENSG\d*)"', line)
    if ref_gene_id is not None:
        line = re.sub(r'MSTRG\.\d*', ref_gene_id.groups(1)[0], line)
        replaced_lines.append(line)
    else:
        replaced_lines.append(line)

gtf_with_repalcement.writelines(replaced_lines)