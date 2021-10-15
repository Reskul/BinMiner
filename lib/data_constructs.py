from lib import Sequence
import numpy as np


class MarkerGene:
    def __init__(self, name):
        self.MG_name = name
        self.contigs = []

    def __str__(self):
        return f"{self.MG_name}:{self.contigs}"

    def __repr__(self):
        return f"{self.MG_name}:{self.contigs}"

    def add_contig(self, contig):
        self.contigs.append(contig)

    def add_contigs(self, contigs):
        for c in contigs:
            self.contigs.append(c)

    def __contains__(self, item) -> bool:
        for c in self.contigs:
            if c == item:
                return True
        return False


class Contig:
    def __init__(self, seq: Sequence, organism=None):
        self.CONTIG_name = seq.header.contig
        self.mgs = []
        self.organism = organism
        self.coverage = None
        self.sequence = seq
        self.kmere_counts = None

    def __str__(self):
        return f"{self.CONTIG_name}/{self.organism}[{self.coverage}]:{self.mgs}"

    def __repr__(self):
        return f"{self.CONTIG_name}/{self.organism}[{self.coverage}]:{self.mgs}"

    def add_mg(self, marker_gene):
        self.mgs.append(marker_gene)

    def __contains__(self, item) -> bool:
        for mg in self.mgs:
            if mg == item:
                return True
        return True

    def to_fasta(self):
        # cut sequences in 80 char pieces or squeeze in \n at pos 80
        length = len(self.sequence)
        n_lines = int(np.floor(length / 80))
        overhang = length % 80
        seq = self.sequence
        seq_lines = []
        for i in range(n_lines):
            seq_lines.append(seq[(i-1)*80:i*80])
        seq_lines.append(seq[length - overhang:])
        seq = '\n'.join(seq_lines)
        # print(n_lines, overhang)
        # print(length, len(seq))
        return f">{self.CONTIG_name};organism={self.organism};cov={self.coverage};mg={self.mgs}{seq}"
