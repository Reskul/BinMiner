class FastaHeaderUniProtKB:
    associated_ma_seq = None

    def __init__(self, db, u_id, entry, prot, org, o_id, gene, ex, ver):
        self.database = db
        self.unique_id = u_id
        self.entry_name = entry
        self.protein_name = prot
        self.organism = org
        self.organism_id = o_id
        self.gene = gene
        self.protein_existence = ex
        self.sequence_version = ver

    def __str__(self):
        return str(self.unique_id)

    def __repr__(self):
        return str(self.unique_id)

    def get_info(self):
        print(self.database, self.unique_id, self.entry_name, self.protein_name, self.organism, self.organism_id,
              self.gene, self.protein_existence, self.sequence_version)

    def add_associated_ma_seq(self, seq_id):
        self.associated_ma_seq = seq_id


class FastaHeaderNCBI:
    def __init__(self, identifier, prot_name, host):
        self.unique_id = identifier
        self.prot_name = prot_name
        self.organism = host

    def __str__(self):
        return "" + str(self.unique_id) + " | " + str(self.prot_name) + " | " + str(self.organism)

    def __repr__(self):
        return "" + str(self.unique_id) + " | " + str(self.prot_name) + " | " + str(self.organism)


class FastaHeaderPRODIGAL:
    def __init__(self, contig, left, right, strand, id, partial, start_type, rbs_motif, rbs_spacer, gc_cont):
        self.contig_part = contig
        self.left = left
        self.right = right
        self.strand = strand
        self.id = id
        self.partial = partial,
        self.start_type = start_type
        self.rbs_motif = rbs_motif
        self.rbs_spacer = rbs_spacer
        self.gc_cont = gc_cont

        self.contig_pure = contig.split('_')[0]

    def __str__(self):
        return f"{self.contig_pure} ID={self.id}"

    def __repr__(self):
        return f"{self.contig_pure} ID={self.id}"


class FastaHeaderMYCC:
    def __init__(self, contig):
        self.contig = contig

    def __str__(self):
        return f"{self.contig}"

    def __repr__(self):
        return f"{self.contig}"
