import ntpath
import os
import subprocess
import numpy as np

from kpal.klib import Profile

from PyQt5.QtCore import *
from sklearn.manifold import TSNE
from lib import FastaReader, Contig, MarkerGene

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class DataLoadingRunnable(QRunnable):
    def __init__(self, *args, **kwargs):
        # called_by, contig_path, kmere_path, datadir, fetchmg_respath, debug=False, create=True
        # create: False -> use fetchmg results
        # create: True -> calc fetchMG
        super().__init__()
        self.call = args[0]
        self.contig_path = args[1]
        self.coverage_path = args[2]
        self.DATADIR = args[3]
        self.perplexity = args[4]
        self.plotstate = args[5]
        self.kmere_path = kwargs.get('kmere_path', None)
        self.DEBUG = kwargs.get('debug', False)
        self.testdata_path = kwargs.get('testdata_path', None)
        self.create = kwargs.get('create', False)
        self.corr = 0
        if self.testdata_path is None:
            self.TEST_MODE = False
        else:
            self.TEST_MODE = True

        dataset_name = ntpath.basename(self.contig_path)

        if not self.create:
            self.fetchmg_respath = args[6]
        else:
            self.fetchmg_binpath = args[6]
            self.prodigal_binpath = args[7]

            self.protein_filename = f"{dataset_name}_prot.fasta"
            self.mg_output_dir = f"mgs_{dataset_name}"

    def run(self) -> None:
        if self.create:
            if not os.path.exists(f"{self.DATADIR}/prodigal/"):
                print(f"Creating {self.DATADIR}/prodigal/")
                os.makedirs(f"{self.DATADIR}/prodigal/")

            completed_prodigal = subprocess.run(
                [self.prodigal_binpath, "-a", f"{self.DATADIR}/prodigal/{self.protein_filename}", "-i",
                 self.contig_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            completed_prodigal.check_returncode()

            completed_fetchmg = subprocess.run(
                [f"./{self.fetchmg_binpath}", "-d", self.contig_path, "-o", f"{self.DATADIR}/{self.mg_output_dir}",
                 "-m extraction", f"{self.DATADIR}/prodigal/{self.protein_filename}"], stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            # completed_fetchmg.check_returncode()  # fetchMG doesn't like being run like this

        kmere_counts_debug = None
        if self.DEBUG and False: #need a new flag, maybe include in test flag but thats not the intention, so deactivating this with false
            # not used atm
            kmere_counts_debug = self.read_kmer_data()

        contigs, mgs = self.read_data()

        kmere_counts = self.calculate_kmere_data()

        failed = False
        plotdata = None

        if self.check_correlation(kmere_counts, contigs):
            # zip Kmere_counts into Contig
            i = 0
            while i < len(contigs):
                contigs[i].kmere_counts = kmere_counts[i]
                i += 1

            # add single dimension coverage representation to contig
            covs_multidim = [c.coverage for c in contigs]
            coverage_dim = len(contigs[0].coverage)
            coverages_1d = None
            if coverage_dim > 1:
                scaled_cov = StandardScaler().fit_transform(covs_multidim)
                coverages_1d = PCA(n_components=1, random_state=5).fit_transform(scaled_cov)
            i = 0
            while i < len(contigs):
                if coverage_dim > 1:
                    contigs[i].coverage_1d = coverages_1d[i]
                else:
                    contigs[i].coverage_1d = contigs[i].coverage
                i += 1

            if self.plotstate == 0:  # Plot Kmere Data
                kmere_sums = np.sum(kmere_counts, 1)
                x_mat = (kmere_counts.T / kmere_sums).T  # Norm data (relative frequencies)
                # T-SNE Data Dim reduction
                plotdata = TSNE(n_components=2, perplexity=self.perplexity).fit_transform(x_mat)
            elif self.plotstate == 1:  # Plot Coverage Data
                # Only if dimension is 2 oder higher
                if len(contigs[0].coverage) > 1:
                    plotdata = TSNE(n_components=2, perplexity=self.perplexity).fit_transform(covs_multidim)
                else:
                    failed = True
                    QMetaObject.invokeMethod(self.call, "data_failed", Qt.QueuedConnection,
                                             Q_ARG(str, "Coverage is 1-Dimensional. Use different Method"))

            elif self.plotstate == 2:  # Plot Combined Data
                n = len(contigs)
                cov_dim = len(contigs[0].coverage)
                kmere_dim = len(contigs[0].kmere_counts)
                x_mat = np.empty((n, kmere_dim + cov_dim))

                for i in range(n):
                    cov = contigs[i].coverage
                    kmere = contigs[i].kmere_counts
                    row = np.hstack((kmere, cov))
                    x_mat[i] = row
                # TODO normalize x_mat here to not over- or underweight coverage
                x_mat_sums = np.sum(x_mat, 1)
                x_mat = (x_mat.T / x_mat_sums).T
                plotdata = TSNE(n_components=2, perplexity=self.perplexity).fit_transform(x_mat)
            else:
                failed = True
                QMetaObject.invokeMethod(self.call, "data_failed", Qt.QueuedConnection, Q_ARG(str, "Unknown Plotstate"))

        else:
            failed = True
            QMetaObject.invokeMethod(self.call, "data_failed", Qt.QueuedConnection,
                                     Q_ARG(str, f"Correlation failed ({self.corr})"))

        if not failed:
            QMetaObject.invokeMethod(self.call, "data_ready", Qt.QueuedConnection, Q_ARG(type(contigs), contigs),
                                     Q_ARG(type(mgs), mgs), Q_ARG(type(plotdata), plotdata))

    def read_data(self):
        test_data = None
        if self.TEST_MODE:
            # read Test Data [contig_name,heritage]
            test_file = open(self.testdata_path, 'r')
            test_data = self.read_test_data(test_file)

        if self.create:
            mg_path = os.path.join(self.DATADIR, self.mg_output_dir)
        else:
            mg_path = self.fetchmg_respath
        cog_files = [file for file in os.listdir(mg_path) if
                     os.path.isfile(os.path.join(mg_path, file)) and file.endswith(".faa")]
        if self.DEBUG:
            print(f"Anzahl: {len(cog_files)}\n{cog_files}")

        # READ-IN Marker Genes from fetchMG results ----------
        # TODO: re-think markergene reading, there is a table file which maybe makes everything easier
        reader = FastaReader(FastaReader.PRODIGAL)
        mgs = []
        for f in cog_files:
            path = os.path.join(mg_path, f)
            header = reader.read_header_only(open(path, 'r'))
            contigs = [c.contig_pure for c in header]
            mg = MarkerGene(f.split('.')[0])
            mg.add_contigs(contigs)
            mgs.append(mg)

        # READ-IN Data from Fasta to get all existing Contigs ----------
        reader = FastaReader(FastaReader.MYCC)
        sequences = reader.read_full_fasta(open(self.contig_path, 'r'))

        coverage_file = open(self.coverage_path, 'r')
        coverage_tup = self.read_coverage(coverage_file)

        contigs = np.empty(len(sequences), dtype=Contig)
        i_idx = 0
        # TODO is there a better method? faster? --> maybe with the table of markergenes but only maybe
        for s in sequences:
            if self.TEST_MODE:
                organism = 'ERROR'
                for i in range(len(test_data)):
                    if s.header.contig == test_data[i][0]:
                        organism = test_data[i][1]
                        break
                c = Contig(s, organism)  # May break sometimes? well.. if files are wrong
            else:
                c = Contig(s)
            for mg in mgs:
                if mg.__contains__(c.CONTIG_name):
                    c.add_mg(mg.MG_name)
            for entry in coverage_tup:
                if c.CONTIG_name == entry[0]:
                    c.coverage = np.array([float(val) for val in entry[1:]], dtype=float)
            if self.TEST_MODE:
                for i in range(len(test_data)):
                    if c.CONTIG_name == test_data[i][0]:
                        c.organism = test_data[i][1]
            contigs[i_idx] = c
            i_idx += 1

        if self.DEBUG:
            print(f"[DEBUG] DataLoadingRunnable.read_dna_data(): Contig Example: {contigs[0]}\t{type(contigs[0])}")
            print(f"[DEBUG] DataLoadingRunnable.read_dna_data(): Number fo Contigs: {len(contigs)} | Type: {type(contigs)} ")
            # print(f"[DEBUG] DataLoadingRunnable.read_dna_data(): Fasta Representation:\n{contigs[0].to_fasta()}")

        mgs = np.array(mgs, dtype=MarkerGene)

        return contigs, mgs

    def calculate_kmere_data(self):
        if self.DEBUG:
            print("[DEBUG] DataLoadingRunnable.run(): Calculating Kmere counts")

        # copied code from P.Meinicke
        # %% kmer profiles
        profs = Profile.from_fasta_by_record(open(self.contig_path), 4)
        p_list = [p for p in profs]
        n = len(p_list)
        p_mat = np.zeros((n, 256))

        for i in range(n):
            p = p_list[i]
            p.balance()
            p_mat[i] = p.counts

        i_vec = -np.ones(256)
        for i in range(256):
            j = p.reverse_complement(i)
            if i_vec[j] == -1:
                i_vec[i] = i

        inds = i_vec[i_vec > -1]
        x_raw_mat = p_mat[:, i_vec > -1]

        return x_raw_mat

    def read_kmer_data(self):
        if self.DEBUG:
            print("[DEBUG] DataLoadingRunnable.run()", self.kmere_path)
        data_raw = np.load(self.kmere_path)
        # dir_path = ntpath.dirname(self.kmere_path)
        # dataset_name = ntpath.basename(self.kmere_path).split('_')[0]
        # access_name = "_".join((dataset_name, "access.npy"))
        # access_path = os.path.join(dir_path, access_name)
        # access_labels = np.load(access_path)
        if self.DEBUG:
            # print(f"[DEBUG] Access Object:{access_labels[0]}\tType:{type(access_labels[0])}")
            print(f"[DEBUG] Datenpunkte:{np.shape(data_raw)}")
            # print(f"[DEBUG] Access Object:{access_labels[0]}")
        return data_raw

    def read_coverage(self, cov_file):
        lines = cov_file.read().split('\n')
        lines_nbr = len(lines)
        lines = lines[:lines_nbr - 1]

        cov_dim = len(lines[0].split('\t'))  # TODO: check if this -1 is needed or not
        data_len = len(lines)
        if self.DEBUG:
            print(f"[DEBUG] DataLoadingRunnable.read_coverage(): CovDim:{cov_dim}")

        covs = np.zeros((1, cov_dim), dtype=str)
        for l in lines:
            res = l.split('\t')
            covs = np.vstack((covs, np.array(res, dtype=str)))  # [res[0], [float(val) for val in res[1:]]]

        covs = covs[1:]
        if self.DEBUG:
            print(f"[DEBUG] DataLoadingRunnable.read_coverage(): Cov's:{covs[0]} Type:{type(covs)}")

        return covs

    def read_test_data(self, test_file):
        entries = test_file.read().split('\n')
        if self.DEBUG:
            tmp = entries[0].split('\t')
            print(f"[DEBUG]DataLoadingRunnable.read_test_data(): {entries[0]}|{tmp}|{len(entries)}")
        arr = np.empty((len(entries), 2), dtype=object)
        contig_name = []
        org_name = []
        for i in range(len(entries)):
            if entries[i].__contains__('\t'):
                e = entries[i].split('\t')
                # contig_name.append(e[0].strip())
                # org_name.append(e[1].strip())
                arr[i, 0] = e[0]
                arr[i, 1] = e[1]
        if self.DEBUG:
            print(f"[DEBUG]DataLoadingRunnable.read_test_data(): {arr.shape}|{arr[0][0]},{arr[0][1]}")
            # print(f"[DEBUG]DataLoadingRunnable.read_test_data(): {len(contig_name)},{len(org_name)}|{contig_name[0]},{org_name[0]}")

        return arr

    def check_correlation(self, kmeres, contigs) -> bool:
        if len(kmeres) == len(contigs):
            n = len(kmeres)
            i = 0
            c_lengths = np.empty(n)
            k_lengths = np.empty(n)
            while i < n:
                contig = contigs[i]
                kmere = kmeres[i]
                c_lengths[i] = len(contig.sequence)
                k_lengths[i] = sum(kmere)
                i += 1
            c_sigma = np.std(c_lengths)
            c_mean = np.mean(c_lengths)
            k_sigma = np.std(k_lengths)
            k_mean = np.mean(k_lengths)
            c_k_cov = 1 / n * sum((c_lengths - c_mean) * (k_lengths - k_mean))
            self.corr = c_k_cov / (c_sigma * k_sigma)
            if self.DEBUG:
                print(
                    f"[DEBUG] DataLoadingRunnable.check_correlation(): Correlation:{self.corr}")
            if self.corr > 0.9:
                return True
            else:
                return False
        else:
            print(f"[ERROR] DataLoadingRunnable.check_correlation(): Kmere and Contig Counts do not match.")
            return False
