class MarkerGene:
    MG_name = None
    contigs = []

    def __init__(self, name):
        self.MG_name = name

    def add_contig(self, contig):
        self.contigs.append(contig)

    def __str__(self):
        return f"{self.MG_name}:{self.contigs}"

    def __repr__(self):
        return f"{self.MG_name}:{self.contigs}"


class Contig:
    CONTIG_name = None
    mgs = []

    def __init__(self, name):
        self.CONTIG_name = name

    def add_mg(self, marker_gene):
        self.mgs.append(marker_gene)

    def __str__(self):
        return f"{self.CONTIG_name}:{self.mgs}"

    def __repr__(self):
        return f"{self.CONTIG_name}:{self.mgs}"
