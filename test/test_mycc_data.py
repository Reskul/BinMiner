from lib import *

reader = FastaReader(FastaReader.MYCC)

contigs = "/home/rom/Dokumente/BA_Data/mycc/10s/10s.fasta"
depths = "/home/rom/Dokumente/BA_Data/mycc/10s/10s.depth.txt"

depths_file = open(depths, 'r')


contigs_file = open(contigs, 'r')
header = reader.read_raw_file(contigs_file)

print(header)
