import ntpath
import os
import subprocess

import numpy as np
from PyQt5.QtCore import *
from sklearn.manifold import TSNE

from fasta import *


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
        contig_lengths = np.sum(data_raw, 1)
        x_mat = (data_raw.T / contig_lengths).T  # Norm data ??? ASK!
        data = TSNE(n_components=2).fit_transform(x_mat)  # T-SNE Data Dim reduction
        print(np.shape(data))
        QMetaObject.invokeMethod(self.call, "set_data", Qt.QueuedConnection, Q_ARG(type(data), data))


class FastaLoadingRunnable(QRunnable):

    def __init__(self, fasta_path, called_by, prodigal_path, fetchmg_path, datadir, debug: bool = False,
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

        contigs = self.read_mgs()
        QMetaObject.invokeMethod(self.call, "protdata_ready", Qt.QueuedConnection, Q_ARG(type(contigs), contigs))

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
            print(f"[DEBUG] Contig Example: {contigs[0]}")
            print(f"[DEBUG] Number fo Contigs: {len(contigs)} | Type: {type(contigs)} ")

        return contigs
