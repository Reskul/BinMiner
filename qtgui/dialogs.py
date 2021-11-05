import ntpath
import numpy as np
import os
from datetime import datetime

from PyQt5.QtWidgets import *
from lib import *

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

matplotlib.rc('font', size=25)


class BinInfoDialog(QDialog):
    STURGES = 0
    SCOTT = 1
    FREEDMAN_DIACONIS = 2

    def __init__(self, parent, debug=False):
        super().__init__(parent)
        self.setWindowTitle("Selected Bin")
        self.DEBUG = debug
        self.parent = parent

        self.selected_kmer = None
        self.selected_cov = None

        self.count_arr = np.array([])
        self.cut_quartiles = False

        self.cut_ckbox = None

        self.figure = None
        self.c_ax = None
        self.k_ax = None

        self.init_gui()

    def init_gui(self):
        # Maybe insert some buttons for output selected sequences as file and stuff

        # histogramm groupbox
        histo_gbox = QGroupBox("Histograms of Coverage and K-mer PC1")
        self.figure = Figure(figsize=(20, 10), dpi=50)
        self.c_ax = self.figure.add_subplot(121)
        self.k_ax = self.figure.add_subplot(122)
        canvas = FigureCanvas(self.figure)
        save_btn = QPushButton("save figure as")
        cut_btn = QPushButton("Cut Quartiles")
        self.cut_ckbox = QCheckBox("Cut Quartiles")

        save_btn.clicked.connect(self.save_clicked)
        cut_btn.clicked.connect(self.cut_clicked)
        self.cut_ckbox.clicked.connect(self.cut_clicked)

        grid = QGridLayout()
        grid.addWidget(canvas, 0, 0, 1, 6)
        grid.addWidget(save_btn, 1, 5)
        grid.addWidget(self.cut_ckbox, 1, 0)
        histo_gbox.setLayout(grid)

        hbox = QHBoxLayout()
        hbox.addWidget(histo_gbox)

        self.setLayout(hbox)
        self.update_gui()

    def save_clicked(self):
        # TODO open file selector
        now = datetime.now().strftime("%d-%m-%y_%H-%M-%S")
        file_name = f"cov_histo{now}.png"
        filepath, _ = QFileDialog.getSaveFileName(self, "Save File", file_name)
        if filepath:
            self.figure.savefig(fname=filepath, dpi=100)

    def cut_clicked(self):
        if self.cut_ckbox.isChecked():
            self.cut_quartiles = True
        else:
            self.cut_quartiles = False
        self.update_gui()

    # def update_selected(self, selected: np.ndarray):
    #     self.selected = selected
    #     if self.DEBUG:
    #         print(f"[DEBUG] BinInfoDialog.update_selected()")
    #
    #     if self.isVisible():
    #         self.update_gui()
    #
    # def set_contigs(self, contigs: Contig):
    #     self.contigs = contigs
    #     if self.DEBUG:
    #         print(f"[DEBUG] BinInfoDialog.set_contigs()\n"
    #               f"\tContigs:{contigs}")
    #
    # def set_markergenes(self, mgs):
    #     if self.DEBUG:
    #         print(f"[DEBUG] Len MGS:{len(mgs)}")
    #     # form into dictionary
    #     for i_idx in range(len(mgs)):
    #         self.mg_dict[mgs[i_idx].MG_name] = i_idx
    #
    #     if self.DEBUG:
    #         print(f"[DEBUG] BinInfoDialog.set_markergenes()\n"
    #               f"\tMG_Dictionary:{self.mg_dict}")
    #     return mgs
    #
    # def set_contigs_and_markergenes(self, contigs, mgs):
    #     self.contigs = contigs
    #     return contigs, self.set_markergenes(mgs)
    #
    # # def find_selected_contigs(self):
    # #     if self.selected is not None:
    # #         contigs = self.contigs[self.selected]
    # #
    # #         if self.DEBUG:
    # #             print(f"[DEBUG] BinInfoDialog.find_selected_contigs()\n"
    # #                   f" Selected Contigs[0]:{contigs[0]}")
    # #
    # #         self.sel_contigs = contigs
    #
    # def calc_values(self):
    #     sel_contigs = self.contigs[self.selected]
    #     coverages = []
    #     kmer_counts = []
    #     self.count_arr = np.zeros(len(self.mgs), dtype=int)  # counting which markergenes exist
    #     for c in sel_contigs:
    #         c_mgs = c.mgs
    #         for mg in c_mgs:
    #             self.count_arr[self.mg_dict[mg]] += 1
    #         coverages.append(c.coverage)
    #         kmer_counts.append(c.kmere_counts)
    #
    #     val_greater_zero = [val > 0 for val in self.count_arr]
    #     completeness = sum(val_greater_zero) / len(self.mgs)
    #     if self.DEBUG:
    #         print(f"[DEBUG] BinInfoDialog.calc_values()\n"
    #               f"\tCounted MG's:{self.count_arr}\n"
    #               f"\tValues greater than 0:{val_greater_zero}\n"
    #               f"\tCompleteness:{completeness}")
    #     max_cnt = max(self.count_arr)
    #     i_max = 2
    #     contamination = 0
    #     while i_max <= max_cnt:
    #         existing = sum([val == i_max for val in self.count_arr])
    #         contamination += existing * (i_max - 1)
    #         i_max += 1
    #     contamination = contamination / len(self.mgs)
    #
    #     if self.DEBUG:
    #         print(f"[DEBUG] BinInfoDialog.calc_values(): Contamination:{contamination}")
    #
    #     n_components = 1
    #
    #     kmer_counts = np.array(kmer_counts, dtype=int)
    #     scaled_kmere_cnt = StandardScaler().fit_transform(kmer_counts)
    #     kmer_counts = PCA(n_components=n_components, random_state=5).fit_transform(scaled_kmere_cnt)
    #
    #     coverages = np.array(coverages, dtype=float)
    #     if len(coverages[0]) > 1:
    #         n_components = 1
    #         scaled_cov = StandardScaler().fit_transform(coverages)
    #         coverages = PCA(n_components=n_components, random_state=5).fit_transform(coverages)
    #     sorted_cov = np.sort(coverages)
    #     if self.cut_quartiles:
    #         q1_idx = int(np.round(len(sorted_cov) * 0.25))
    #         q3_idx = int(np.round(len(sorted_cov) * 0.75))
    #         return completeness, contamination, sorted_cov[q1_idx:q3_idx], kmer_counts[q1_idx:q3_idx]
    #     else:
    #         return completeness, contamination, sorted_cov, kmer_counts

    def update_data(self, sel_kmer_counts, sel_coverages):
        n_components = 1
        n_selected = len(sel_coverages)

        sel_kmer_counts = np.array(sel_kmer_counts, dtype=int)
        sel_coverages = np.array(sel_coverages, dtype=float)

        if n_selected == 1:
            sel_kmer_counts.reshape(1, -1)
            sel_coverages.reshape(1, -1)

        scaled_kmere_cnt = StandardScaler().fit_transform(sel_kmer_counts)
        sel_kmer_counts = PCA(n_components=n_components, random_state=5).fit_transform(scaled_kmere_cnt)

        if len(sel_coverages[0]) > 1:
            n_components = 1
            scaled_cov = StandardScaler().fit_transform(sel_coverages)
            sel_coverages = PCA(n_components=n_components, random_state=5).fit_transform(sel_coverages)

        self.selected_cov = np.sort(sel_coverages)
        self.selected_kmer = np.sort(sel_kmer_counts)

        self.update_gui()

    def update_histo(self, cov=None, kmer=None, bin_rule=0):
        if cov is None and kmer is None:
            cov = self.selected_cov
            kmer = self.selected_kmer

        n = len(cov)
        cov_sigma = np.std(cov)
        kmere_sigma = np.std(cov)
        q3, q1 = np.percentile(cov, [75, 25])
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.update_histo: Quartiles:{q1, q3}; Sigma:{cov_sigma}")
        if bin_rule == 0:
            # Sturges-Rule for Bin number
            k_bins = c_bins = 1 + np.log2(n)
        elif bin_rule == 1:
            # Scott-Rule
            c_bins = (3.49 * cov_sigma) / np.cbrt(n)  # std:StandardAbweichung Sigma | cbrt:cubicroot
            k_bins = (3.49 * kmere_sigma) / np.cbrt(n)
        elif bin_rule == 2:
            # Freedman-Diaconis Rule
            k_bins = c_bins = (2 * (q3 - q1)) / np.cbrt(n)
        else:
            k_bins = c_bins = n / 5

        c_bins_round = int(np.round(c_bins))
        k_bins_round = int(np.round(k_bins))

        if self.DEBUG:
            print(
                f"[DEBUG] BinInfoDialog.update_histo(): N:{n}; BinsRaw:{c_bins}; Bins:{c_bins_round}; Rule:{bin_rule}")

        # if c_bins_round < 1:
        #     c_bins = int(np.round(n / 5))
        #     print(f"[DEBUG] BinInfoDialog.update_histo(): Bins:{c_bins}")

        self.c_ax.cla()
        self.k_ax.cla()

        self.c_ax.hist(cov, bins=c_bins_round)
        self.k_ax.hist(kmer, bins=k_bins_round)

        self.c_ax.set_title("Coverage of selected contigs")
        self.k_ax.set_title("K-mer distribution of selected contigs")

        self.c_ax.set_xlabel("Coverage")
        self.c_ax.set_ylabel("Frequency")
        self.k_ax.set_xlabel("K-mer PC1")
        self.k_ax.set_ylabel("Frequency")

        self.figure.canvas.draw_idle()

    def show(self) -> None:
        if self.selected_cov is not None and self.selected_kmer is not None:
            super().show()
            self.update_gui()
        elif self.DEBUG:
            print(f"[ERROR] Something is missing. You need to insert Contigs ans Markergenes first.")

    def update_gui(self):
        if self.isVisible() and self.selected_cov is not None and self.selected_kmer is not None:
            if self.cut_quartiles:
                n_selected = len(self.selected_cov)
                q1_idx = int(np.round(n_selected * 0.25))
                q3_idx = int(np.round(n_selected * 0.75))
                self.update_histo(cov=self.selected_cov[q1_idx:q3_idx], kmer=self.selected_kmer[q1_idx:q3_idx])
            else:
                self.update_histo()


class NameSelectedDialog(QDialog):

    def __init__(self, parent, debug=False):
        super().__init__(parent)
        self.setWindowTitle("Name Organism")
        self.setModal(True)
        self.DEBUG = debug
        self.parent = parent

        self.name_input: QLineEdit = None
        self.ok_btn = None
        self.cancel_btn = None

        self.init_gui()

    def init_gui(self):
        title_text = QLabel("Name the prototype genome")
        self.name_input = QLineEdit()
        self.ok_btn = QPushButton("OK")
        self.ok_btn.clicked.connect(self.ok_clicked)
        self.cancel_btn = QPushButton("CANCEL")
        self.cancel_btn.clicked.connect(self.cancel_clicked)

        layout = QGridLayout()
        layout.addWidget(title_text, 0, 1)
        layout.addWidget(self.name_input, 1, 1)
        layout.addWidget(self.cancel_btn, 2, 0, 1, 1)
        layout.addWidget(self.ok_btn, 2, 2, 1, 1)

        self.setLayout(layout)

    def ok_clicked(self):
        name = self.name_input.text()
        name.strip()
        if name != '':
            self.parent.save_to_file(name)
        else:
            self.parent.save_to_file()

        self.close()

    def cancel_clicked(self):
        self.close()
