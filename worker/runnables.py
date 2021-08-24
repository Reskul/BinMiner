import ntpath
import os
import subprocess

import numpy as np
from PyQt5.QtCore import *
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from lib import *


class LoadingNpyRunnable(QRunnable):
    def __init__(self, data_path, called_by, debug=False):
        super().__init__()
        self.path = data_path
        self.call = called_by
        self.DEBUG = debug

    def run(self):
        if self.DEBUG:
            print("[DEBUG] LoadingNpyRunnable.run()", self.path)
        data_raw = np.load(self.path)
        dir_path = ntpath.dirname(self.path)
        dataset_name = ntpath.basename(self.path).split('_')[0]
        access_name = "_".join((dataset_name, "access.npy"))
        access_path = os.path.join(dir_path, access_name)
        access_labels = np.load(access_path)

        contig_lengths = np.sum(data_raw, 1)
        x_mat = (data_raw.T / contig_lengths).T  # Norm data ??? ASK!
        data = TSNE(n_components=2).fit_transform(x_mat)  # T-SNE Data Dim reduction
        if self.DEBUG:
            print(f"[DEBUG] Access Object:{access_labels[0]}\tType:{type(access_labels[0])}")
            print(f"[DEBUG] Datenpunkte:{np.shape(data)}")
            print(f"[DEBUG] Access Object:{access_labels[0]}")
        QMetaObject.invokeMethod(self.call, "set_data", Qt.QueuedConnection, Q_ARG(type(data), data),
                                 Q_ARG(type(access_labels), access_labels))


class FastaLoadingRunnable(QRunnable):

    def __init__(self, called_by, fasta_path, datadir, prodigal_path=None, fetchmg_path=None, debug: bool = False,
                 only_analyze: bool = False, contig_translation: bool = False, translation_file=None):
        super().__init__()
        self.path = fasta_path
        self.call = called_by
        self.prodigal = prodigal_path
        self.fetchMG = fetchmg_path
        self.DATADIR = datadir
        if contig_translation:
            self.translation_path = translation_file

        self.DEBUG = debug
        self.ONLY_ANALYZE = only_analyze
        self.TRANSLATE_CONTIGS = contig_translation

        fasta_filename = ntpath.basename(fasta_path)
        dataset_name, filetype = fasta_filename.split('.', 1)
        self.protein_file = f"{dataset_name}_prot.fasta"
        self.mg_output_dir = f"mgs_{dataset_name}"

    def run(self) -> None:
        if self.DEBUG:
            print(
                f"[DEBUG] FastaLoadingRunnable.run()\nOnly-Analyze:{self.ONLY_ANALYZE}\n{self.path}\n{self.prodigal}\n{self.fetchMG}")

        if not self.ONLY_ANALYZE:

            if not os.path.exists(f"{self.DATADIR}/prodigal/"):
                print(f"Creating {self.DATADIR}/prodigal/")
                os.makedirs(f"{self.DATADIR}/prodigal/")

            completed_prodigal = subprocess.run(
                [self.prodigal, "-a", f"{self.DATADIR}/prodigal/{self.protein_file}", "-i", self.path],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            completed_prodigal.check_returncode()

            completed_fetchmg = subprocess.run(
                ["perl", self.fetchMG, "-d", self.path, "-o", f"{self.DATADIR}/{self.mg_output_dir}", "-m extraction",
                 f"{self.DATADIR}/prodigal/{self.protein_file}"],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            # completed_fetchmg.check_returncode()

        contigs, mgs = self.read_mgs()
        QMetaObject.invokeMethod(self.call, "protdata_ready", Qt.QueuedConnection, Q_ARG(type(contigs), contigs),
                                 Q_ARG(type(mgs), mgs))

    def read_mgs(self):
        mg_path = os.path.join(self.DATADIR, self.mg_output_dir)
        fasta_path = self.path
        cog_files = [file for file in os.listdir(mg_path) if
                     os.path.isfile(os.path.join(mg_path, file)) and file.endswith(".faa")]
        if self.DEBUG:
            print(f"Anzahl: {len(cog_files)}\n{cog_files}")
        # READ-IN Marker Genes from fetchMG results
        reader = FastaReader(FastaReader.PRODIGAL)
        mgs = []
        for f in cog_files:
            path = os.path.join(mg_path, f)
            header = reader.read_raw_file(open(path, 'r'))
            contigs = [c.contig_pure for c in header]
            mg = MarkerGene(f.split('.')[0])
            mg.add_contigs(contigs)
            mgs.append(mg)

        # if self.DEBUG:
        #     print(f"[DEBUG] MG[0]:\n{mgs[0]}")

        # READ-IN Data from Fasta to get all existing Contigs
        reader = FastaReader(FastaReader.MYCC)
        header = reader.read_raw_file(open(fasta_path, 'r'))
        contigs = np.empty(len(header), dtype=Contig)
        content = None
        if self.TRANSLATE_CONTIGS:
            trans_file = open(self.translation_path, 'r')
            content = trans_file.read()
        i_idx = 0
        # TODO is there a better method? faster?
        for h in header:
            c = Contig(h.contig)
            for mg in mgs:
                if mg.__contains__(h.contig):
                    c.add_mg(mg.MG_name)
            if self.TRANSLATE_CONTIGS:
                x = content.find(c.CONTIG_name)
                y = content.find("\t", x)
                z = content.find("\n", y)
                y += 1
                c.REAL_name = content[y:z]

            contigs[i_idx] = c
            i_idx += 1

        if self.DEBUG:
            print(f"[DEBUG] Contig Example: {contigs[0]}\t{type(contigs[0])}")
            print(f"[DEBUG] Number fo Contigs: {len(contigs)} | Type: {type(contigs)} ")

        mgs = np.array(mgs)

        return contigs, mgs


class DataLoadingRunnable(QRunnable):
    def __init__(self, *args, **kwargs):
        # called_by, contig_path, kmere_path, datadir, fetchmg_respath, debug=False, create=True
        # create: False -> use fetchmg results
        # create: True -> calc fetchMG
        super().__init__()
        self.call = args[0]
        self.contig_path = args[1]
        self.coverage_path = args[2]
        self.kmere_path = args[3]
        self.DATADIR = args[4]
        self.DEBUG = kwargs.get('debug', False)
        self.create = kwargs.get('create', False)

        dataset_name = ntpath.basename(self.contig_path)

        if not self.create:
            self.fetchmg_respath = args[5]
        else:
            self.fetchmg_binpath = args[5]
            self.prodigal_binpath = args[6]

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
                ["perl", self.fetchmg_binpath, "-d", self.contig_path, "-o", f"{self.DATADIR}/{self.mg_output_dir}",
                 "-m extraction", f"{self.DATADIR}/prodigal/{self.protein_filename}"], stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            # completed_fetchmg.check_returncode()  # fetchMG doesn't like being run like this

        datapoints = self.read_data()
        contigs, mgs = self.read_mgs()
        QMetaObject.invokeMethod(self.call, "data_ready", Qt.QueuedConnection, Q_ARG(type(contigs), contigs),
                                 Q_ARG(type(mgs), mgs), Q_ARG(type(datapoints), datapoints))

    def read_mgs(self):
        if self.create:
            mg_path = os.path.join(self.DATADIR, self.mg_output_dir)
        else:
            mg_path = self.fetchmg_respath
        cog_files = [file for file in os.listdir(mg_path) if
                     os.path.isfile(os.path.join(mg_path, file)) and file.endswith(".faa")]
        if self.DEBUG:
            print(f"Anzahl: {len(cog_files)}\n{cog_files}")

        # READ-IN Marker Genes from fetchMG results ----------
        reader = FastaReader(FastaReader.PRODIGAL)
        mgs = []
        for f in cog_files:
            path = os.path.join(mg_path, f)
            header = reader.read_raw_file(open(path, 'r'))
            contigs = [c.contig_pure for c in header]
            mg = MarkerGene(f.split('.')[0])
            mg.add_contigs(contigs)
            mgs.append(mg)

        # READ-IN Data from Fasta to get all existing Contigs ----------
        reader = FastaReader(FastaReader.MYCC)
        header = reader.read_raw_file(open(self.contig_path, 'r'))

        coverage_file = open(self.coverage_path, 'r')
        coverage_tup = self.read_coverage(coverage_file)

        contigs = np.empty(len(header), dtype=Contig)
        i_idx = 0
        # TODO is there a better method? faster?
        for h in header:
            c = Contig(h.contig)
            for mg in mgs:
                if mg.__contains__(h.contig):
                    c.add_mg(mg.MG_name)
            for entry in coverage_tup:  # TODO Again, too slow. better way?
                if c.CONTIG_name == entry[0]:
                    c.coverage = entry[1]
            contigs[i_idx] = c
            i_idx += 1

        if self.DEBUG:
            print(f"[DEBUG] Contig Example: {contigs[0]}\t{type(contigs[0])}")
            print(f"[DEBUG] Number fo Contigs: {len(contigs)} | Type: {type(contigs)} ")

        mgs = np.array(mgs)

        return contigs, mgs

    def read_data(self):
        if self.DEBUG:
            print("[DEBUG] DataLoadingRunnable.run()", self.kmere_path)
        data_raw = np.load(self.kmere_path)
        # dir_path = ntpath.dirname(self.kmere_path)
        # dataset_name = ntpath.basename(self.kmere_path).split('_')[0]
        # access_name = "_".join((dataset_name, "access.npy"))
        # access_path = os.path.join(dir_path, access_name)
        # access_labels = np.load(access_path)
        # TODO: Contig Längen mit den Counts abgleichen
        contig_lengths = np.sum(data_raw, 1)
        x_mat = (data_raw.T / contig_lengths).T  # Norm data ??? ASK!
        data = TSNE(n_components=2).fit_transform(x_mat)  # T-SNE Data Dim reduction
        if self.DEBUG:
            # print(f"[DEBUG] Access Object:{access_labels[0]}\tType:{type(access_labels[0])}")
            print(f"[DEBUG] Datenpunkte:{np.shape(data)}")
            # print(f"[DEBUG] Access Object:{access_labels[0]}")
        return data

    def read_coverage(self, cov_file):
        lines = cov_file.read().split('\n')
        lines_nbr = len(lines)
        lines = lines[:lines_nbr - 1]

        cov_dim = len(lines[0].split('\t')) - 1
        data_len = len(lines)

        names = []
        covs = np.empty((1, cov_dim))
        for l in lines:
            res = l.split('\t')
            tup = res[0]
            names.append(tup)
            covs = np.vstack((covs, [float(val) for val in res[1:]]))

        # PCA on multidim Coverages
        n_components = 1  # wont be changed because Histogramm can only show 1-dim
        covs = covs[1:]
        if cov_dim > 1:
            if self.DEBUG:
                print(f"[DEBUG] DataLoadingRunnable.read_coverage(): Applying PCA")
            # Standard Scaler
            scaled_data = StandardScaler().fit_transform(covs)  # Maybe use this if non-scaled fails
            reduced_covs = PCA(n_components, random_state=5).fit_transform(covs)
        else:
            reduced_covs = covs

        i = 0
        collection = []
        while i < data_len:
            collection.append((names[i], reduced_covs[i]))
            i += 1
        return np.array(collection, dtype=tuple)
