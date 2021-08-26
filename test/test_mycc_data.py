from lib import *

reader = FastaReader(FastaReader.MYCC)

contigs = "/home/rom/Dokumente/BA_Data/mycc/10s/10s.fasta"
depths = "/home/rom/Dokumente/BA_Data/mycc/10s/10s.depth.txt"

assembly_25s = "/home/rom/Dokumente/BA_Data/mycc/25s/assembly.fa"

# depths_file = open(depths, 'r')


# contigs_file = open(contigs, 'r')
assembly_file = open(assembly_25s,'r')
header = reader.read_header_only(assembly_file)
test = reader.read_full_fasta(assembly_file)
print(header)
