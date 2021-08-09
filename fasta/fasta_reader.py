import numpy as np
from fasta.fasta_header import *
from fasta.sequence import Sequence


class FastaReader:
    UNIPROTKB = 1
    NCBI = 2

    def __init__(self, db):
        self.db = db  # db: Database 1=UniProtKB, 2=NCBI
        pass

    def read_raw_file(self, collection_file):
        data_vec = collection_file.read().split('>')[1:]
        header = np.empty(len(data_vec), dtype=FastaHeaderUniProtKB)
        if self.db == self.UNIPROTKB:
            i_idx = 0
            for entry in data_vec:
                info_line = entry.split("\n")[0]
                info_line = info_line.split("|")
                db = info_line[0]
                unique_id = info_line[1]
                info_line = info_line[2]

                # identifier indices/marker
                org_name_idx = info_line.find("OS=")
                org_id_idx = info_line.find("OX=")
                gene_name_idx = info_line.find("GN=")
                prot_existence_idx = info_line.find("PE=")
                seq_ver_idx = info_line.find("SV=")

                # information
                prefix = info_line[:org_name_idx].split()
                entry_name = prefix[0]
                prot_name = " ".join(prefix[1:])
                org_name = info_line[org_name_idx + 3:org_id_idx].strip()
                org_id = info_line[org_id_idx + 3:gene_name_idx].strip()
                gene_name = info_line[gene_name_idx + 3:prot_existence_idx].strip()
                prot_existence = info_line[prot_existence_idx + 3:seq_ver_idx].strip()
                seq_ver = info_line[seq_ver_idx + 3:].strip()

                header[i_idx] = FastaHeaderUniProtKB(db, unique_id, entry_name, prot_name, org_name, org_id, gene_name, prot_existence, seq_ver)
                i_idx = i_idx + 1
        else:
            pass
        return header

    # read in multiple alignment as fasta
    def read_ma_fasta(self, ma_file):
        # MA processing
        data = ma_file.read()
        data = data.split(">")[1:]
        seq_cnt = len(data)
        seq_vec = np.empty(seq_cnt, dtype=FastaHeaderUniProtKB)
        if self.db == self.UNIPROTKB:
            i_idx = 0
            for item in data:
                seq_fasta = item.split("\n")

                name = seq_fasta[:1][0]

                seq = seq_fasta[1:]
                seq = "".join(seq)
                seq_vec[i_idx] = Sequence(seq, name)
                i_idx = i_idx + 1
        else:
            pass
        return seq_vec

    # compare raw seq header with ma data
    def compare_raw_ma(self, raw_header_vec, ma_seq_vec):
        if len(raw_header_vec) != len(ma_seq_vec):
            print("[INPUT_CHECK]Length Error")
            return False
        else:
            ma_seq_names_vec = [seq.unique_id for seq in ma_seq_vec]
            for header in raw_header_vec:
                if header.unique_id not in ma_seq_names_vec:
                    print("[INPUT_CHECK]Name Error")
                    return False
        return True

    def combine_sequence_data(self, header_vec, sequence_vec):
        """Combines the information of both sequence files to gather all information in one place"""
        # check if raw sequence data fits to the MA data
        if not self.compare_raw_ma(header_vec, sequence_vec):
            print("Sequence File and Multiple Alignment file do not fit.")
            return None

        # add header to sequence object
        only_uid_header_vec = np.array([h.unique_id for h in header_vec])

        for seq in sequence_vec:
            old_idx = np.where(only_uid_header_vec == seq.unique_id)
            seq.add_header(header_vec[old_idx])
        return sequence_vec

    def read(self, ma_filename, collection_filename):
        header_vec = self.read_raw_file(open(collection_filename))
        ma_vec = self.read_ma_fasta(open(ma_filename))
        seq_vec = self.combine_sequence_data(header_vec, ma_vec)

        return seq_vec
