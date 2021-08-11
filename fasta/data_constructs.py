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
    def __init__(self, name, real_name=None):
        self.CONTIG_name = name
        self.mgs = []
        self.REAL_name = real_name

    def __str__(self):
        return f"{self.CONTIG_name}/{self.REAL_name}:{self.mgs}"

    def __repr__(self):
        return f"{self.CONTIG_name}/{self.REAL_name}:{self.mgs}"

    def add_mg(self, marker_gene):
        self.mgs.append(marker_gene)

    def __contains__(self, item) -> bool:
        for mg in self.mgs:
            if mg == item:
                return True
        return True
