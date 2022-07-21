from lib import FastaReader

if __name__ == '__main__':
    contig_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Ergebnisse\\190622_10s\\Bin00_n2185.fasta", 'r')

    reader = FastaReader(FastaReader.BINMINER)
    sequences = reader.read_full_fasta(contig_file)
    yes = 0
    no = 0
    for seq in sequences:
        if len(seq) < 1000:
            no += 1
            print(seq.header.organism)
        else:
            yes += 1
    print(yes, no)

    count = 0
    for seq in sequences:
        if seq.header.organism == 'ERROR':
            print(seq.header.contig, len(seq))
            count += 1
    print('Error count:', count)
