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

    def __init__(self, parent, selected: np.ndarray = None, contigs=None, mgs=None, debug=False):
        super().__init__(parent)
        self.setWindowTitle("Selected Bin")
        self.DEBUG = debug
        self.parent = parent
        self.mg_dict = {}
        self.selected = selected
        self.contigs = contigs
        self.count_arr = np.array([])
        self.sel_contigs = None
        self.mgs = self.set_markergenes(mgs)

        self.cut_quartiles = False

        self.completeness_nbr_lbl = None
        self.containment_nbr_lbl = None
        self.cut_ckbox = None

        self.figure = None
        self.c_ax = None
        self.k_ax = None

        self.init_gui()

    def init_gui(self):
        numbers_gbox = QGroupBox("Statistics")
        completeness_lbl = QLabel("Completeness:")
        self.completeness_nbr_lbl = QLabel()

        containment_lbl = QLabel("Contamination:")
        self.containment_nbr_lbl = QLabel()

        form_layout = QFormLayout()
        form_layout.addRow(completeness_lbl, self.completeness_nbr_lbl)
        form_layout.addRow(containment_lbl, self.containment_nbr_lbl)

        numbers_gbox.setLayout(form_layout)

        # histogramm groupbox
        histo_gbox = QGroupBox("Coverage Histogram")
        self.figure = Figure(figsize=(20, 10), dpi=50)
        self.c_ax = self.figure.add_subplot(121)
        self.k_ax = self.figure.add_subplot(122)
        canvas = FigureCanvas(self.figure)
        save_btn = QPushButton("Save As")
        cut_btn = QPushButton("Cut Quartiles")
        self.cut_ckbox = QCheckBox("Cut Quartiles")

        save_btn.clicked.connect(self.save_clicked)
        cut_btn.clicked.connect(self.cut_clicked)
        self.cut_ckbox.clicked.connect(self.cut_clicked)

        grid = QGridLayout()
        grid.addWidget(canvas, 0, 0, 1, 3)
        grid.addWidget(save_btn, 1, 2)
        grid.addWidget(self.cut_ckbox, 1, 0)
        histo_gbox.setLayout(grid)

        hbox = QHBoxLayout()
        hbox.addWidget(numbers_gbox)
        hbox.addWidget(histo_gbox)

        self.setLayout(hbox)
        self.update_gui()

    def save_clicked(self):
        now = datetime.now().strftime("%d-%m-%y_%H-%M-%S")
        file_name = f"cov_histo{now}"
        d_path = os.path.join(self.parent.cfg.homepath, 'coverage_histograms')
        if not os.path.exists(d_path):
            os.makedirs(d_path)
        f_path = os.path.join(self.parent.cfg.homepath, 'coverage_histograms', file_name)
        self.figure.savefig(fname=f_path, dpi=100)

    def cut_clicked(self):
        if self.cut_ckbox.isChecked():
            self.cut_quartiles = True
        else:
            self.cut_quartiles = False
        self.update_gui()

    def update_selected(self, selected: np.ndarray):
        self.selected = selected
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.update_selected()")

        if self.isVisible():
            self.update_gui()

    def set_contigs(self, contigs: Contig):
        self.contigs = contigs
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.set_contigs()\n"
                  f"\tContigs:{contigs}")

    def set_markergenes(self, mgs):
        if self.DEBUG:
            print(f"[DEBUG] Len MGS:{len(mgs)}")
        # form into dictionary
        for i_idx in range(len(mgs)):
            self.mg_dict[mgs[i_idx].MG_name] = i_idx

        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.set_markergenes()\n"
                  f"\tMG_Dictionary:{self.mg_dict}")
        return mgs

    def set_contigs_and_markergenes(self, contigs, mgs):
        self.contigs = contigs
        return contigs, self.set_markergenes(mgs)

    # def find_selected_contigs(self):
    #     if self.selected is not None:
    #         contigs = self.contigs[self.selected]
    #
    #         if self.DEBUG:
    #             print(f"[DEBUG] BinInfoDialog.find_selected_contigs()\n"
    #                   f" Selected Contigs[0]:{contigs[0]}")
    #
    #         self.sel_contigs = contigs

    def calc_values(self):
        # print(f"DIGGAH {self.sel_contigs}")
        sel_contigs = self.contigs[self.selected]
        coverages = []
        kmere_counts = []
        self.count_arr = np.zeros(len(self.mgs), dtype=int)
        for c in sel_contigs:
            c_mgs = c.mgs
            # print(f"LOL {c}")
            for mg in c_mgs:
                # print(f"ROFL {mg}")
                self.count_arr[self.mg_dict[mg]] += 1
            coverages.append(c.coverage)
            kmere_counts.append(c.kmere_counts)

        val_greater_zero = [val > 0 for val in self.count_arr]
        completeness = sum(val_greater_zero) / len(self.mgs)
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.calc_values()\n"
                  f"\tCounted MG's:{self.count_arr}\n"
                  f"\tValues greater than 0:{val_greater_zero}\n"
                  f"\tCompleteness:{completeness}")
        max_cnt = max(self.count_arr)
        i_max = 2
        contamination = 0
        while i_max <= max_cnt:
            existing = sum([val == i_max for val in self.count_arr])
            contamination += existing * (i_max - 1)
            i_max += 1
        contamination = contamination / len(self.mgs)

        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.calc_values(): Contamination:{contamination}")

        n_components = 1

        kmere_counts = np.array(kmere_counts, dtype=int)
        scaled_kmere_cnt = StandardScaler().fit_transform(kmere_counts)
        kmere_counts = PCA(n_components=n_components, random_state=5).fit_transform(kmere_counts)

        coverages = np.array(coverages, dtype=float)
        if len(coverages[0]) > 1:
            n_components = 1
            scaled_cov = StandardScaler().fit_transform(coverages)
            coverages = PCA(n_components=n_components, random_state=5).fit_transform(coverages)
        sorted_cov = np.sort(coverages)
        if self.cut_quartiles:
            q1_idx = int(np.round(len(sorted_cov) * 0.25))
            q3_idx = int(np.round(len(sorted_cov) * 0.75))
            return completeness, contamination, sorted_cov[q1_idx:q3_idx], kmere_counts[q1_idx:q3_idx]
        else:
            return completeness, contamination, sorted_cov, kmere_counts

    def update_histo(self, cov, kmere, bin_rule=0):
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
        self.k_ax.hist(kmere, bins=k_bins_round)

        self.c_ax.set_title("Coverage of selected contigs")
        self.k_ax.set_title("K-mere distribution of selected contigs")

        self.c_ax.set_xlabel("Coverage")
        self.c_ax.set_ylabel("Frequency")
        self.k_ax.set_xlabel("K-mere Count")
        self.k_ax.set_ylabel("Frequency")

        self.figure.canvas.draw_idle()

    def show(self) -> None:
        if self.contigs is not None and self.mgs is not None:
            super().show()
            self.update_gui()
        elif self.DEBUG:
            print(f"[ERROR] Something is missing. You need to insert Contigs ans Markergenes first.")

    def update_gui(self):
        if self.selected is not None and self.isVisible() and self.contigs is not None and self.mgs is not None:
            completeness, contamination, coverages, kmere_counts = self.calc_values()

            self.update_histo(coverages, kmere_counts, self.STURGES)
            self.completeness_nbr_lbl.setText(str(completeness))
            self.containment_nbr_lbl.setText(str(contamination))
