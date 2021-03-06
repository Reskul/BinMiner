import numpy as np
from lib.fasta_header import *
from lib.sequence import Sequence


class FastaReader:
    UNIPROTKB = 1
    NCBI = 2
    PRODIGAL = 3
    MYCC = 4
    BINMINER = 5
    RAYMETA = 6

    def __init__(self, origin):
        self.origin = origin  # origin/Database 1=UniProtKB, 2=NCBI, 3=Prodigal

    def read_header_only(self, collection_file):
        data_vec = collection_file.read().split('>')[1:]
        header_vec = None
        if self.origin == self.UNIPROTKB:
            header_vec = np.empty(len(data_vec), dtype=FastaHeaderUniProtKB)
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

                header_vec[i_idx] = FastaHeaderUniProtKB(db, unique_id, entry_name, prot_name, org_name, org_id, gene_name,
                                                         prot_existence, seq_ver)
                i_idx = i_idx + 1
        elif self.origin == self.PRODIGAL:
            header_vec = np.empty(len(data_vec), dtype=FastaHeaderPRODIGAL)
            i_idx = 0
            for entry in data_vec:
                info_line = entry.split('\n')[0]
                first = info_line.split('#')
                second = first[4]
                first = first[:4]
                second = second.split(';')

                name = first[0].strip()
                left = first[1].strip()
                right = first[2].strip()
                strand = first[3].strip()
                id = second[0].split('=')[1].strip()
                partial = second[1].split('=')[1]
                start_type = second[2].split('=')[1]
                rbs_motif = second[3].split('=')[1]
                rbs_spacer = second[4].split('=')[1]
                gc_cont = second[5].split('=')[1]

                header_vec[i_idx] = FastaHeaderPRODIGAL(name, left, right, strand, id, partial, start_type, rbs_motif,
                                                        rbs_spacer, gc_cont)
                i_idx += 1
        elif self.origin == self.MYCC:
            header_vec = np.empty(len(data_vec), dtype=FastaHeaderMYCC)
            i_idx = 0
            for entry in data_vec:
                info_line = entry.split('\n')[0]
                header_vec[i_idx] = FastaHeaderMYCC(info_line.strip())
                i_idx += 1
        elif self.origin == self.BINMINER:
            header_vec = np.empty(len(data_vec), dtype=FastaHeaderBINMINER)
            i_idx = 0
            while i_idx < len(data_vec):
                info_line = data_vec[i_idx].split('\n')[0]
                info_vec = info_line.split(';')
                contig = info_vec[0].strip()
                organism = info_vec[1].strip().split('=')[1]
                coverage = info_vec[2].strip()
                mgs = info_vec[3].strip()

                header_vec[i_idx] = FastaHeaderBINMINER(contig, organism, coverage, mgs)
                i_idx += 1
        else:
            print(f"[ERROR] Unknown header. {self.origin}")
        return header_vec

    # read in multiple alignment as fasta
    def read_full_fasta(self, fasta_file):
        # MA processing
        data_vec = fasta_file.read()
        data_vec = data_vec.split(">")[1:]
        seq_cnt = len(data_vec)
        seq_vec = None
        if self.origin == self.UNIPROTKB:
            seq_vec = np.empty(seq_cnt, dtype=FastaHeaderUniProtKB)
            i_idx = 0
            for item in data_vec:
                seq_fasta = item.split("\n")

                name = seq_fasta[:1][0]

                seq = seq_fasta[1:]
                seq = "".join(seq)
                seq_vec[i_idx] = Sequence(seq, name)
                i_idx = i_idx + 1
        elif self.origin == self.MYCC:
            seq_vec = np.empty(seq_cnt, dtype=Sequence)
            i_idx = 0
            for entry in data_vec:
                parts = entry.split('\n')
                info_line = parts[0]
                seq = "".join(parts[1:])
                h = FastaHeaderMYCC(info_line.split(' ')[0])

                seq_vec[i_idx] = Sequence(seq, i_idx, header=h)
                i_idx += 1
        elif self.origin == self.RAYMETA:
            seq_vec = np.empty(seq_cnt, dtype=Sequence)
            i_idx = 0
            for entry in data_vec:
                parts = entry.split('\n')
                info_line = parts[0]
                seq = "".join(parts[1:])
                h = FastaHeaderRAYMETA(info_line.split(' ')[0], info_line.split(' ')[1])

                seq_vec[i_idx] = Sequence(seq, i_idx, header=h)
                i_idx += 1
        elif self.origin == self.BINMINER:
            seq_vec = np.empty(seq_cnt, dtype=Sequence)
            i_idx = 0
            for entry in data_vec:
                parts = entry.split('\n')
                info_line = parts[0]
                seq = "".join(parts[1:])
                info_vec = info_line.split(';')
                name = info_vec[0]
                org = info_vec[1].split('=')[1]
                cov = [float(x) for x in info_vec[2].split('=')[1].strip('[]').split(',')]
                mgs = info_vec[3].split('=')[1].strip('[]').split(',')
                h = FastaHeaderBINMINER(name, org, cov, mgs)

                seq_vec[i_idx] = Sequence(seq, i_idx, header=h)
                i_idx += 1
        else:
            print(f"[ERROR] FastaReader.read_full_fasta(): Unknown Header.")
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

    def read_ma_and_header(self, ma_filename, collection_filename):
        header_vec = self.read_header_only(open(collection_filename))
        ma_vec = self.read_full_fasta(open(ma_filename))
        seq_vec = self.combine_sequence_data(header_vec, ma_vec)

        return seq_vec


if __name__ == '__main__':
    file = open('/home/rom/Documents/MGB/prodigal/10s_prot.fasta')
    print(type(file))
    reader = FastaReader(FastaReader.PRODIGAL)
    header = reader.read_header_only(file)
    print(header[0])
