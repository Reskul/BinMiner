from lib import FastaReader

if __name__ == '__main__':
    contig_in_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_contigs.fasta")
    labels_in_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_labels.txt")
    cov_in_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s.depth.txt")

    new_contigs_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_cutoff.fasta", "w")
    removed_contigs_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_contigs_removed.fasta", "w")

    reader = FastaReader(FastaReader.MYCC)
    sequences = reader.read_full_fasta(contig_in_file)
    contig_in_file.close()

    cutoff = 1000
    new_contigs = []
    removed_contigs = []
    for seq in sequences:
        if len(seq) < cutoff:
            removed_contigs.append(seq.header.contig)
            print(seq.header.contig)
            removed_contigs_file.write(f">{seq.header.contig}\n")
            removed_contigs_file.write(f"{seq}\n")
            removed_contigs_file.flush()
        else:
            new_contigs.append(seq)
            new_contigs_file.write(f">{seq.header.contig}\n")
            new_contigs_file.write(f"{seq}\n")
            new_contigs_file.flush()

    new_contigs_file.close()
    removed_contigs_file.close()

    print(f"Removed:{len(removed_contigs)} | Remaining:{len(new_contigs)}")

    lines = labels_in_file.readlines()
    cov_lines = cov_in_file.readlines()
    labels_in_file.close()
    new_lines = []
    new_cov = []
    i = 0
    while i < len(lines):
        contig, organism = lines[i].split('\t')
        contig_cov, cov = cov_lines[i].split('\t')
        if not removed_contigs.__contains__(contig):
            new_lines.append(lines[i])
        if not removed_contigs.__contains__(contig_cov):
            new_cov.append(cov_lines[i])
        i += 1

    new_cov_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_cutoff.depth.txt", "w")
    new_cov_file.writelines(new_cov)
    new_cov_file.flush()
    new_cov_file.close()

    new_labels_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\10s_cutoff_labels.txt", "w")
    new_labels_file.writelines(new_lines)
    new_labels_file.flush()
    new_labels_file.close()
