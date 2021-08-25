import ntpath
import numpy as np
import os
from datetime import date
from datetime import datetime

from PyQt5.QtWidgets import *
from lib import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class PathsDialog(QDialog):
    OS = None
    PARENT = None
    DEFAULTPATH = None
    DEBUG = None

    cfg = None

    prodigal_path = None
    fetchMG_path = None
    datadir_path = None

    def __init__(self, parent: QMainWindow, defaultpath, cfg: Configurator, os: str = 'Linux', debug=False) -> None:
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle("Set Paths")
        self.setFixedWidth(900)
        self.cfg = cfg
        # set Variables ----------
        self.OS = os
        self.PARENT = parent
        self.DEFAULTPATH = defaultpath
        self.DEBUG = debug

        # Entries ----------
        prodigal_lbl = QLabel("Prodigal:")
        prodigal_btn = QPushButton("Wählen")
        prodigal_btn.clicked.connect(self.set_prodigal)
        self.prodigal_le = QLineEdit()
        self.prodigal_le.setReadOnly(True)  # TODO: Necessary?

        fetchmg_lbl = QLabel("FetchMG:")
        fetchmg_btn = QPushButton("Wählen")
        fetchmg_btn.clicked.connect(self.set_fetchmg)
        self.fetchMG_le = QLineEdit()
        self.fetchMG_le.setReadOnly(True)  # TODO: Necessary?

        data_lbl = QLabel("Daten:")
        data_btn = QPushButton("Wählen")
        data_btn.clicked.connect(self.set_datadir)
        self.data_le = QLineEdit()
        self.data_le.setReadOnly(True)  # TODO: Necessary?

        # TEMPORARY PATHS ----------
        mg_lbl = QLabel("FetchMG Ergebnisse:")
        mg_btn = QPushButton("Wählen")
        mg_btn.clicked.connect(self.set_fetchmg_results)
        self.mg_le = QLineEdit()
        self.mg_le.setReadOnly(True)

        fasta_lbl = QLabel("Original Contig Daten:")
        fasta_btn = QPushButton("Wählen")
        fasta_btn.clicked.connect(self.set_fasta_data)
        self.fasta_le = QLineEdit()
        self.fasta_le.setReadOnly(True)

        fasta_names_lbl = QLabel("Contig Übersetzung(Optional):")
        fasta_names_btn = QPushButton("Wählen")
        fasta_names_btn.clicked.connect(self.set_fasta_translation)
        self.fasta_names_le = QLineEdit()
        self.fasta_names_le.setReadOnly(True)

        start_btn = QPushButton("Berechnen")
        start_btn.clicked.connect(self.start_result_processing)

        # Control btns
        ok_btn = QPushButton("Ok")
        ok_btn.clicked.connect(self.ok_clicked)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.cancel_clicked)

        # Set Texts if path existing ----------
        self.prodigal_le.setText(self.cfg.read(self.cfg.PRODIGAL_KEY))
        self.fetchMG_le.setText(self.cfg.read(self.cfg.FETCHMG_KEY))
        self.data_le.setText(self.cfg.read(self.cfg.DATA_KEY))

        # Design Separator ----------
        line = QFrame()
        line.setGeometry(0, 0, 200, 3)
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)

        # Overall Layout ----------
        layout = QGridLayout(self)
        layout.addWidget(prodigal_lbl, 0, 0)
        layout.addWidget(self.prodigal_le, 0, 1)
        layout.addWidget(prodigal_btn, 0, 2)
        layout.addWidget(fetchmg_lbl, 1, 0)
        layout.addWidget(self.fetchMG_le, 1, 1)
        layout.addWidget(fetchmg_btn, 1, 2)
        # layout.addWidget(data_lbl, 2, 0)
        # layout.addWidget(self.data_le, 2, 1)
        # layout.addWidget(data_btn, 2, 2)

        layout.addWidget(line, 2, 0, 1, 3)
        layout.addWidget(QLabel("Temporär"), 3, 0)

        layout.addWidget(mg_lbl, 4, 0)
        layout.addWidget(self.mg_le, 4, 1)
        layout.addWidget(mg_btn, 4, 2)

        layout.addWidget(fasta_lbl, 5, 0)
        layout.addWidget(self.fasta_le, 5, 1)
        layout.addWidget(fasta_btn, 5, 2)

        layout.addWidget(fasta_names_lbl, 6, 0)
        layout.addWidget(self.fasta_names_le, 6, 1)
        layout.addWidget(fasta_names_btn, 6, 2)

        layout.addWidget(start_btn, 7, 2)

        layout.addWidget(ok_btn, 8, 2)
        layout.addWidget(cancel_btn, 8, 0)

    def set_prodigal(self):
        path = None
        if self.OS == 'Windows':
            path, _ = QFileDialog.getOpenFileName(self, 'Select Prodigal Executable', self.DEFAULTPATH,
                                                  'Executables (*.exe)')
        elif self.OS == 'Linux':
            path, _ = QFileDialog.getOpenFileName(self, 'Select Prodigal Binary', self.DEFAULTPATH,
                                                  'Binaries (*.linux)')
        # TODO add OSX solution somehow
        elif self.DEBUG:
            print("[DEBUG] Window.set_prodigal(): WHY U USING OSX?")
        if path:
            if self.DEBUG:
                print(f"[DEBUG] Window.set_prodigal(): path = {path}")
            self.prodigal_le.setText(path)
            self.PARENT.prodigal_path = path
            self.prodigal_path = path

    def set_fetchmg(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select FetchMG Script', self.DEFAULTPATH,
                                              'Perl Scripts (*.pl)')
        if path:
            self.fetchMG_le.setText(path)
            self.PARENT.fetchMG_path = path
            self.fetchMG_path = path

    def set_datadir(self):
        """select dir where all data is stored"""
        dir_path = None
        if self.OS == 'Linux':
            dir_path = QFileDialog.getExistingDirectory(self, 'Select Directory', '/home', QFileDialog.ShowDirsOnly)
        elif self.OS == 'Windows':
            dir_path = QFileDialog.getExistingDirectory(self, 'Select Directory', 'C:\\Users', QFileDialog.ShowDirsOnly)
        if dir_path:
            self.data_le.setText(dir_path)
            self.PARENT.DATADIR = dir_path
            self.datadir_path = dir_path

    def set_fetchmg_results(self):
        path = QFileDialog.getExistingDirectory(self, 'Select Directory', self.cfg.homepath, QFileDialog.ShowDirsOnly)
        if path:
            self.mg_le.setText(path)
            # self.PARENT.TMP_fetchMG_results_path = path

    def set_fasta_data(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Fasta File', self.cfg.homepath, "Fasta Files (*.fasta)")
        if path:
            self.fasta_le.setText(path)
            # self.PARENT.TMP_prodigal_results_path = path

    def set_fasta_translation(self):
        if ntpath.exists(self.fasta_le.text()):
            in_path = ntpath.dirname(self.fasta_le.text())
        else:
            in_path = self.cfg.homepath

        path, _ = QFileDialog.getOpenFileName(self, 'Select Translation', in_path, "SPE Text Files (*.spe.txt)")

        if path:
            self.fasta_names_le.setText(path)

    def start_result_processing(self):
        self.PARENT.process_markergenes(self.mg_le.text(), self.fasta_le.text(), self.fasta_names_le.text())

    def ok_clicked(self):
        if self.datadir_path:
            self.cfg.write(self.cfg.DATA_KEY, self.datadir_path)
        if self.prodigal_path:
            self.cfg.write(self.cfg.PRODIGAL_KEY, self.prodigal_path)
        if self.fetchMG_path:
            self.cfg.write(self.cfg.FETCHMG_KEY, self.fetchMG_path)
        self.close()

    def cancel_clicked(self):
        self.close()


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

        self.figure = None
        self.ax = None

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
        self.figure = Figure(figsize=(16, 9), dpi=45)
        self.ax = self.figure.add_subplot(111)
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
        self.figure.savefig(fname=f_path)

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

        if self.contigs is not None and self.mgs is not None:

            self.find_selected_contigs()

        if self.isVisible():
            self.update_gui()

    def set_contigs(self, contigs: Contig):
        self.contigs = contigs
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.set_contigs()\n"
                  f"\tContigs:{contigs}")

        if self.selected is not None and self.mgs is not None:
            self.find_selected_contigs()

    def set_markergenes(self, mgs):
        if self.DEBUG:
            print(f"[DEBUG] Len MGS:{len(mgs)}")
        # form into dictionary
        for i_idx in range(len(mgs)):
            self.mg_dict[mgs[i_idx].MG_name] = i_idx

        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.set_markergenes()\n"
                  f"\tMG_Dictionary:{self.mg_dict}")
        if self.selected is not None and self.contigs is not None:
            self.find_selected_contigs()
        return mgs

    def set_contigs_and_markergenes(self, contigs, mgs):
        self.contigs = contigs
        return contigs, self.set_markergenes(mgs)

    def find_selected_contigs(self):
        if self.selected is not None:
            contigs = self.contigs[self.selected]

            if self.DEBUG:
                print(f"[DEBUG] BinInfoDialog.find_selected_contigs()\n"
                      f" Selected Contigs[0]:{contigs[0]}")

            self.sel_contigs = contigs

    def calc_values(self):
        # print(f"DIGGAH {self.sel_contigs}")
        coverages = []
        self.count_arr = np.zeros(len(self.mgs), dtype=int)
        for c in self.sel_contigs:
            c_mgs = c.mgs
            # print(f"LOL {c}")
            for mg in c_mgs:
                # print(f"ROFL {mg}")
                self.count_arr[self.mg_dict[mg]] += 1
            coverages.append(c.coverage)

        val_greater_zero = [val > 0 for val in self.count_arr]
        completeness = sum(val_greater_zero) / len(self.mgs)
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.calc_values()\n"
                  f"\tCounted MG's:{self.count_arr}\n"
                  f"\tValues greater than 0:{val_greater_zero}\n"
                  f"\tCompleteness:{completeness}")
        # completeness = sum([val > 0 for val in self.count_arr])
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

        coverages = np.array(coverages, dtype=float)
        sorted_cov = np.sort(coverages)
        if self.cut_quartiles:
            q1_idx = int(np.round(len(sorted_cov) * 0.25))
            q3_idx = int(np.round(len(sorted_cov) * 0.75))
            print(f'################ {len(sorted_cov[q1_idx:q3_idx])}')
            return completeness, contamination, sorted_cov[q1_idx:q3_idx]
        else:
            return completeness, contamination, sorted_cov

    def update_histo(self, cov, bin_rule=0):
        n = len(cov)
        sigma = np.std(cov)
        q3, q1 = np.percentile(cov, [75, 25])
        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.update_histo: Quartiles:{q1, q3}; Sigma:{sigma}")
        if bin_rule == 0:
            # Sturges-Rule for Bin number
            bins = 1 + np.log2(n)
        elif bin_rule == 1:
            # Scott-Rule
            bins = (3.49 * sigma) / np.cbrt(n)  # std:StandardAbweichung Sigma | cbrt:cubicroot
        elif bin_rule == 2:
            # Freedman-Diaconis Rule
            bins = (2 * (q3 - q1)) / np.cbrt(n)
        else:
            bins = n / 5

        bins_round = int(np.round(bins))

        if self.DEBUG:
            print(f"[DEBUG] BinInfoDialog.update_histo(): N:{n}; BinsRaw:{bins}; Bins:{bins_round}; Rule:{bin_rule}")

        if bins_round < 1:
            bins = int(np.round(n / 5))
            print(f"[DEBUG] BinInfoDialog.update_histo(): Bins:{bins}")

        self.ax.cla()
        # TODO check if cov changes -> should do because values printed out do change
        self.ax.hist(cov, bins=bins_round)
        self.figure.canvas.draw_idle()

    def show(self) -> None:
        if self.contigs is not None and self.mgs is not None:
            super().show()
            self.update_gui()
        elif self.DEBUG:
            print(f"[ERROR] Something is missing. You need to insert Contigs ans Markergenes first.")

    def update_gui(self):
        if self.selected is not None and self.isVisible():
            completeness, contamination, coverages = self.calc_values()

            self.update_histo(coverages, self.STURGES)
            self.completeness_nbr_lbl.setText(str(completeness))
            self.containment_nbr_lbl.setText(str(contamination))
