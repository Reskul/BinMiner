import ntpath

from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from matplotlib import patches
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

from KDEpy import FFTKDE
from skimage.feature import peak_local_max

from lib import *
from .dialogs import BinInfoDialog


class QFileInputLine(QLineEdit):
    clicked = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)

    def is_empty(self):
        if self.text() == "":
            return True
        else:
            return False

    def mousePressEvent(self, a0: QtGui.QMouseEvent) -> None:
        self.clicked.emit()


class InputGUI(QWidget):
    def __init__(self, parent=None, cfg: Configurator = None, debug=False):
        super().__init__(parent)
        # INSTANCE VARIABLES
        self.DEBUG = debug
        self.config = cfg
        self.contig_sequences_path = None
        self.contig_coverage_path = None
        self.kmere_dataset_path = None
        self.fetchmg_result_path = None
        self.prodigal_binpath = None
        self.fetchmg_binpath = None

        # INIT GUI ---------- ---------
        f_layout = QFormLayout()
        form_gbox = QGroupBox()
        general_layout = QVBoxLayout()
        radio_gbox = QGroupBox()
        radio_hbox = QHBoxLayout()

        self.fetchmg_result_form = QFormLayout()
        self.binaries_input_form = QFormLayout()
        self.fetchmg_results_gbox = QGroupBox()
        self.binaries_input_gbox = QGroupBox()

        # Contig DNA Sequence ----------
        self.contig_dataset_le = QFileInputLine()
        self.contig_dataset_le.clicked.connect(self.cd_clicked)
        f_layout.addRow(QLabel("Contig Sequences File"), self.contig_dataset_le)

        # Contig Coverage File
        self.contig_coverage_le = QFileInputLine()
        self.contig_coverage_le.clicked.connect(self.cc_clicked)
        f_layout.addRow((QLabel("Contig Coverage File")), self.contig_coverage_le)

        # K-Mere Data ----------
        self.kmere_dataset_le = QFileInputLine()
        self.kmere_dataset_le.clicked.connect(self.km_clicked)
        f_layout.addRow(QLabel("K-Mere Data File"), self.kmere_dataset_le)

        # TSNE Perplexity Parameter
        perplex_lbl = QLabel("Perplexity of T-SNE")
        self.perplex_spbox = QSpinBox()
        self.perplex_spbox.setMinimum(5)
        self.perplex_spbox.setMaximum(50)
        self.perplex_spbox.setValue(30)

        f_layout.addRow(perplex_lbl, self.perplex_spbox)

        # Radiobutton's Select next Section
        # o FetchMG Ergebnisse | o Selber erstellen
        # >< Ergebnisse Path   | >< Prodigal Path
        #                      | >< FetchMG path
        self.results_radbtn = QRadioButton("FetchMG Results")
        self.source_radbtn = QRadioButton("Calculate MarkerGenes")
        self.results_radbtn.setChecked(True)

        self.results_radbtn.clicked.connect(self.radio_result_clicked)
        self.source_radbtn.clicked.connect(self.radio_source_clicked)

        radio_hbox.addWidget(self.results_radbtn)
        radio_hbox.addWidget(self.source_radbtn)
        radio_hbox.addStretch(1)

        # FetchMG Results Form
        self.fetchmg_res_le = QFileInputLine()
        self.fetchmg_res_le.clicked.connect(self.fmg_clicked)

        self.fetchmg_result_form.addRow(QLabel("FetchMG Results:"), self.fetchmg_res_le)
        self.fetchmg_result_form.addRow(QLabel(""), QLabel(""))

        # Binaries Input Form
        self.prodigal_path_le = QFileInputLine()
        self.fetchmg_path_le = QFileInputLine()
        self.prodigal_path_le.clicked.connect(self.pp_clicked)
        self.fetchmg_path_le.clicked.connect(self.fmgp_clicked)

        self.binaries_input_form.addRow(QLabel("Prodigal Bin"), self.prodigal_path_le)
        self.binaries_input_form.addRow(QLabel("FetchMG Bin"), self.fetchmg_path_le)

        self.binaries_input_gbox.setLayout(self.binaries_input_form)
        self.binaries_input_gbox.setVisible(False)

        # "Next" Button
        next_btn = QPushButton("Next")
        next_btn.clicked.connect(self.next_clicked)

        # Same Procedure as last year Miss Sophie?
        self.last_ckbox = QCheckBox("Use Last Values")
        self.last_ckbox.clicked.connect(self.lastckbox_clicked)

        btn_layout = QGridLayout()
        btn_layout.addWidget(self.last_ckbox, 0, 0, 1, 3)
        btn_layout.addWidget(next_btn, 0, 4, 1, 1)

        # FINALIZE GUI
        form_gbox.setLayout(f_layout)
        radio_gbox.setLayout(radio_hbox)
        self.fetchmg_results_gbox.setLayout(self.fetchmg_result_form)

        general_layout.addStretch(1)
        general_layout.addWidget(form_gbox)
        general_layout.addWidget(radio_gbox)
        general_layout.addWidget(self.fetchmg_results_gbox)
        general_layout.addWidget(self.binaries_input_gbox)
        general_layout.addLayout(btn_layout)
        general_layout.addStretch(1)
        self.setLayout(general_layout)

    def cd_clicked(self):
        # path = QFileDialog.getExistingDirectory(self, 'Select Contig Fasta-Data File', self.config.homepath,
        #                                        QFileDialog.ShowDirsOnly)

        path, _ = QFileDialog.getOpenFileName(self, 'Select Contig Fasta-Data File', self.config.homepath,
                                              'Fasta Files (*.fasta *.faa *.fa *.fna)')

        if path:
            self.contig_dataset_le.setText(path)
            self.contig_sequences_path = path

    def cc_clicked(self):
        if not self.contig_dataset_le.is_empty():
            dirname = ntpath.dirname(self.contig_dataset_le.text())
        else:
            dirname = self.cfg.homepath

        path, _ = QFileDialog.getOpenFileName(self, 'Select Contig Coverage File', dirname,
                                              'Depth/Text Files (*.depth *.depth.txt)')
        if path:
            self.contig_coverage_le.setText(path)
            self.contig_coverage_path = path

    def km_clicked(self):
        # path = QFileDialog.getExistingDirectory(self, 'Select K-mere Data Directory', self.config.homepath,
        #                                         QFileDialog.ShowDirsOnly)
        path, _ = QFileDialog.getOpenFileName(self, 'Select Precomputed K-mere Data', self.config.homepath,
                                              'Numpy Files (*.npy)')
        if path:
            self.kmere_dataset_le.setText(path)
            self.kmere_dataset_path = path

    def radio_result_clicked(self):
        self.radio_changed(True)

    def radio_source_clicked(self):
        self.radio_changed(False)

    def radio_changed(self, flag: bool):
        if flag:
            # Result
            self.binaries_input_gbox.setVisible(False)
            self.fetchmg_results_gbox.setVisible(True)
        else:
            # Source
            self.fetchmg_results_gbox.setVisible(False)
            self.binaries_input_gbox.setVisible(True)

    def fmg_clicked(self):
        path = QFileDialog.getExistingDirectory(self, 'Select FetchMG Result Directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
        if path:
            self.fetchmg_res_le.setText(path)
            self.fetchmg_result_path = path

    def pp_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Prodigal Binary', self.DEFAULTPATH,
                                              'Binaries (*.linux)')
        if path:
            self.prodigal_binpath = path
            self.prodigal_path_le.setText(path)

    def fmgp_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select FetchMG Script', self.DEFAULTPATH,
                                              'Perl Scripts (*.pl)')
        if path:
            self.fetchmg_binpath = path
            self.fetchmg_path_le.setText(path)

    def next_clicked(self):
        if self.contig_sequences_path is not None and self.kmere_dataset_path is not None:
            if self.results_radbtn.isChecked() and self.fetchmg_result_path is not None:
                if self.config:
                    self.config.write(self.config.CONTIGSEQ_KEY, self.contig_sequences_path)
                    self.config.write(self.config.CONTIGCOV_KEY, self.contig_coverage_path)
                    self.config.write(self.config.KMERE_KEY, self.kmere_dataset_path)
                    self.config.write(self.config.FETCHMGRES_KEY, self.fetchmg_result_path)
                self.parent().process_input(self.contig_sequences_path, self.contig_coverage_path,
                                            self.kmere_dataset_path, self.perplex_spbox.value(), fetchmg_respath=self.fetchmg_result_path)

            elif self.source_radbtn.isChecked() and self.fetchmg_binpath is not None and self.prodigal_binpath is not None:
                if self.cfg:
                    self.config.write(self.config.CONTIGSEQ_KEY, self.contig_sequences_path)
                    self.config.write(self.config.CONTIGCOV_KEY, self.contig_coverage_path)
                    self.config.write(self.config.KMERE_KEY, self.kmere_dataset_path)
                    self.config.write(self.config.FETCHMG_KEY, self.fetchmg_binpath)
                    self.config.write(self.config.PRODIGAL_KEY, self.prodigal_binpath)
                self.parent().process_input(self.contig_sequences_path, self.contig_coverage_path,
                                            self.kmere_dataset_path, prodigal_path=self.prodigal_binpath,
                                            fetchmg_path=self.fetchmg_binpath)
            else:
                print(f"[ERROR] Something went terribly wrong with radio buttons.")

    def lastckbox_clicked(self):
        if self.last_ckbox.isChecked():
            # Fill LineEdits with Values
            if not self.config.is_new:
                print(f"[DEBUG] InputGUI.lastckbox_clicked(): Config wurde NICHT neu erstellt.")
                contig_seq = self.config.read(Configurator.CONTIGSEQ_KEY)
                contig_cov = self.config.read(Configurator.CONTIGCOV_KEY)
                kmere_set = self.config.read(Configurator.KMERE_KEY)
                fetchmg_res = self.config.read(Configurator.FETCHMGRES_KEY)
                fetchmg_bin = self.config.read(Configurator.FETCHMG_KEY)
                prodigal_bin = self.config.read(Configurator.PRODIGAL_KEY)

                self.contig_sequences_path = contig_seq
                self.contig_coverage_path = contig_cov
                self.kmere_dataset_path = kmere_set
                self.fetchmg_result_path = fetchmg_res
                self.fetchmg_binpath = fetchmg_bin
                self.prodigal_binpath = prodigal_bin

                self.contig_dataset_le.setText(contig_seq)
                self.contig_coverage_le.setText(contig_cov)
                self.kmere_dataset_le.setText(kmere_set)
                self.fetchmg_res_le.setText(fetchmg_res)
                self.fetchmg_path_le.setText(fetchmg_bin)
                self.prodigal_path_le.setText(prodigal_bin)
        else:
            # clear all LineEdits
            self.contig_dataset_le.setText("")
            self.contig_coverage_le.setText("")
            self.kmere_dataset_le.setText("")
            self.fetchmg_res_le.setText("")
            self.fetchmg_path_le.setText("")
            self.prodigal_path_le.setText("")

            self.contig_sequences_path = None
            self.contig_coverage_path = None
            self.kmere_dataset_path = None
            self.fetchmg_result_path = None
            self.fetchmg_binpath = None
            self.prodigal_binpath = None


class SelectGUI(QWidget):
    def __init__(self, kmere_data, contigs, mgs, parent=None, debug=False):
        super().__init__(parent)
        # VARIABLE Initialization
        self.DEBUG = debug
        self.data = kmere_data
        self.cont_data = None
        self.n_selected = 0
        self.bw = 1
        self.contours = 25
        self.selected_vec = np.zeros(len(self.data))
        self.grid_points = 100

        self.analyze_widget = BinInfoDialog(self.parent(), contigs=contigs, mgs=mgs, debug=self.DEBUG)

        # BACK BUTTON
        back_btn = QPushButton("<- Back")
        back_btn.clicked.connect(self.back_clicked)
        back_btn.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Data Visualization ----------
        fig = Figure()
        self.ax = fig.add_subplot(111)
        self.canvas = FigureCanvas(fig)
        # self.canvas.resize(300, 300)
        toolbar = NavigationToolbar(self.canvas, self)

        diagram_layout = QGridLayout()
        diagram_layout.addWidget(back_btn, 0, 0, 1, 1)
        diagram_layout.addWidget(self.canvas, 1, 0, 1, 3)
        diagram_layout.addWidget(toolbar, 2, 0, 1, 1)

        # Matplotlib Interaction ----------
        self.canvas.mpl_connect("button_press_event", self.on_mpl_press)

        # Process Selected Layout ----------
        process_sel_btn = QPushButton("Ausgewählte Analysieren")
        process_sel_btn.clicked.connect(self.analyze_selected)
        self.sel_lbl = QLabel()
        self.set_selected_cnt()
        selected_layout = QHBoxLayout()
        selected_layout.addWidget(self.sel_lbl)
        selected_layout.addWidget(process_sel_btn)

        diagram_layout.addLayout(selected_layout, 2, 2, 1, 1)

        # Slider Section ----------
        bw_slider_lbl = QLabel("Bandwidth")
        cont_slider_lbl = QLabel("Contours")
        self.bw_nbr_lbl = QLabel(str(self.bw))
        self.cont_nbr_lbl = QLabel(str(self.contours))
        self.bw_slider = QSlider(Qt.Vertical)
        self.cont_slider = QSlider(Qt.Vertical)

        self.bw_slider.setRange(1, 20)
        self.bw_slider.setValue(self.bw)
        self.cont_slider.setRange(1, 50)
        self.cont_slider.setValue(25)

        self.bw_slider.sliderReleased.connect(self.on_bw_slider_change)
        self.cont_slider.sliderReleased.connect(self.on_cont_slider_change)

        slider_layout = QGridLayout()
        slider_layout.addWidget(self.bw_nbr_lbl, 0, 0, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(self.bw_slider, 1, 0, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(bw_slider_lbl, 2, 0)
        slider_layout.addWidget(self.cont_nbr_lbl, 0, 1, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(self.cont_slider, 1, 1, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(cont_slider_lbl, 2, 1)

        # GENERAL LAYOUT
        layout = QHBoxLayout()
        layout.addLayout(diagram_layout)
        layout.addLayout(slider_layout)

        self.update_plot()
        self.setLayout(layout)

    def set_selected_cnt(self):
        if self.selected_vec is not None:
            count = sum(self.selected_vec)
            self.n_selected = count
            self.sel_lbl.setText(f"Ausgewählt: {count}")
        else:
            self.sel_lbl.setText(f"Ausgewählt: --")

    def update_plot(self, highlighted_cont=None, col=None):
        """Updates Matplotlib plots, is called when stuff changed in Data"""
        if self.data is not None:
            self.ax.clear()
            self.ax.tick_params(axis='x', labelsize='14')
            self.ax.tick_params(axis='y', labelsize='14')
            self.ax.patches = []
            x, y, z = self.update_kde()
            if highlighted_cont is not None:
                # TODO: Selection-> Working Visualization: draw only outer line instead of outer&inner of contour level
                patch = patches.PathPatch(highlighted_cont, facecolor=col, lw=1, edgecolor='black', fill=False)
                self.ax.add_patch(patch)
            self.update_datapoints()
            self.update_peaks(x, y, z)
            self.canvas.draw()

    def update_kde(self):
        """Only updates the KDE (density) plot"""
        if self.DEBUG:
            print("[DEBUG] MainWindow.update_kde()")
        grid, points = FFTKDE(norm=1, bw=self.bw).fit(self.data).evaluate(self.grid_points)
        x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
        z = points.reshape(self.grid_points, self.grid_points).T
        # self.cont_data = self.ax.contour(x, y, z, self.contours, cmap='RdBu_r')
        self.cont_data = self.ax.contourf(x, y, z, self.contours, cmap='RdBu_r')
        return x, y, z

    def update_datapoints(self):
        """Only updates the data points"""
        self.ax.scatter(self.data[:, 0], self.data[:, 1], marker=".", s=2, c=self.selected_vec)

    def update_peaks(self, x, y, z):
        """Only updates the peak points"""
        peakdata = peak_local_max(z, threshold_rel=0.02)
        self.ax.scatter(x[peakdata[:, 1]], y[peakdata[:, 0]], marker=10, c="orange")

    def on_bw_slider_change(self):
        self.bw = self.bw_slider.value()
        self.bw_nbr_lbl.setText(str(self.bw))
        self.update_plot()

    def on_cont_slider_change(self):
        self.contours = self.cont_slider.value()
        self.cont_nbr_lbl.setText(str(self.contours))
        self.update_plot()

    def on_mpl_press(self, e):
        """Matplotlib 'press' Event, calculates path and points selected by event"""
        if self.DEBUG:
            print("[DEBUG] MainWindow.on_mpl_press()")
            print("event.xdata", e.xdata)
            print("event.ydata", e.ydata)
        path_list = []
        if self.cont_data is not None:
            for path_collection in self.cont_data.collections:
                check, _ = path_collection.contains(e)
                if check:
                    # path.contains_point()->bool oder path.contains_points()-> bool array
                    paths = path_collection.get_paths()
                    if self.DEBUG:
                        print(len(paths))
                    for p in paths:
                        if p.contains_point((e.xdata, e.ydata)):
                            path_list.append(p)

            last = len(path_list) - 1
            self.selected_vec = np.array(path_list[last].contains_points(self.data), dtype=bool)
            if self.DEBUG:
                print("[DEBUG] SelectGUI.on_mpl_press()")
                print("Nbr of found paths:", len(path_list))
                print(self.selected_vec[0], np.shape(self.selected_vec), '\n', np.shape(self.data))

            self.update_plot(highlighted_cont=path_list[last])  # , col='green')
            self.analyze_widget.update_selected(self.selected_vec)
            self.set_selected_cnt()

    def analyze_selected(self):
        """Takes Selected Datapoints and checks in MG Data for MG's and calculates coverage and contamination"""
        if self.n_selected > 0:
            self.analyze_widget.update_selected(self.selected_vec)

            if not self.analyze_widget.isVisible():
                self.analyze_widget.show()

    def back_clicked(self):
        self.analyze_widget.close()
        self.parent().STATUS = self.parent().STATUS_INPUT
        self.parent().determine_widget()

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        super().closeEvent(a0)
        self.cfg.on_close()
