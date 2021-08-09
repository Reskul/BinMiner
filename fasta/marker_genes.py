class MarkerGene:
    MG_number = None
    contigs = []

    def __init__(self, number):
        self.MG_number = number

    def add_contig(self, contig):
        self.contigs.append(contig)
