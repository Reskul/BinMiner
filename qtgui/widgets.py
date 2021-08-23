from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from matplotlib import patches
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
import os

from KDEpy import FFTKDE
from skimage.feature import peak_local_max

from cfg import *
from lib import *
from .dialogs import BinInfoDialog


class QFileInputLine(QLineEdit):
    clicked = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)

    def mousePressEvent(self, a0: QtGui.QMouseEvent) -> None:
        self.clicked.emit()


class InputGUI(QWidget):
    def __init__(self, parent=None, cfg: Configurator = None, debug=False):
        super().__init__(parent)
        # INSTANCE VARIABLES
        self.DEBUG = debug
        self.config = cfg
        self.contig_dataset_path = None
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
        f_layout.addRow(QLabel("Contig Dataset"), self.contig_dataset_le)

        # K-Mere Data ----------
        self.kmere_dataset_le = QFileInputLine()
        self.kmere_dataset_le.clicked.connect(self.km_clicked)
        f_layout.addRow(QLabel("K-Mere Dataset"), self.kmere_dataset_le)

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

        # FINALIZE GUI
        form_gbox.setLayout(f_layout)
        radio_gbox.setLayout(radio_hbox)
        self.fetchmg_results_gbox.setLayout(self.fetchmg_result_form)

        general_layout.addStretch(1)
        general_layout.addWidget(form_gbox)
        general_layout.addWidget(radio_gbox)
        general_layout.addWidget(self.fetchmg_results_gbox)
        general_layout.addWidget(self.binaries_input_gbox)
        general_layout.addWidget(next_btn, alignment=Qt.AlignRight)
        general_layout.addStretch(1)
        self.setLayout(general_layout)

    def cd_clicked(self):
        # path, _ = QFileDialog.getOpenFileName(self, 'Open Contig File', self.DEFAULTPATH,'Numpy Files (*.fasta)')
        path = QFileDialog.getExistingDirectory(self, 'Select Contig Data Directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
        if path:
            self.contig_dataset_le.setText(path)
            self.contig_dataset_path = path

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
        if self.contig_dataset_path is not None and self.kmere_dataset_path is not None:
            if self.results_radbtn.isChecked() and self.fetchmg_result_path is not None:
                self.parent().process_input(self.contig_dataset_path, self.kmere_dataset_path,
                                            fetchmg_respath=self.fetchmg_result_path)
            elif self.source_radbtn.isChecked() and self.fetchmg_binpath is not None and self.prodigal_binpath is not None:
                self.parent().process_input(self.contig_dataset_path, self.kmere_dataset_path,
                                            prodigal_path=self.prodigal_binpath, fetchmg_path=self.fetchmg_binpath)
            else:
                print(f"[ERROR] Something went terribly wrong with radio buttons.")


class SelectGUI(QWidget):
    def __init__(self, kmere_data, contigs, mgs, parent=None, debug=False):
        super().__init__(parent)
        # VARIABLE Initialization
        self.DEBUG = debug
        self.data = kmere_data
        self.cont_data = None
        self.selected_nbr = 0
        self.bw = 1
        self.contours = 25
        self.selected_data = np.zeros(len(self.data))
        self.grid_points = 100

        self.analyze_widget = BinInfoDialog(self.parent(), contigs=contigs, mgs=mgs, debug=self.DEBUG)

        # Data Visualization ----------
        fig = Figure()
        self.ax = fig.add_subplot(111)
        self.canvas = FigureCanvas(fig)
        # self.canvas.resize(300, 300)
        toolbar = NavigationToolbar(self.canvas, self)

        diagram_layout = QGridLayout()
        diagram_layout.addWidget(self.canvas, 0, 0, 1, 3)
        diagram_layout.addWidget(toolbar, 1, 0, 1, 1)

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

        diagram_layout.addLayout(selected_layout, 1, 2, 1, 1)

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
        if self.selected_data is not None:
            count = sum(self.selected_data)
            self.selected_nbr = count
            self.sel_lbl.setText(f"Ausgewählt: {count}")
        else:
            self.sel_lbl.setText(f"Ausgewählt: --")

    def update_plot(self, highlighted_cont=None, col=None):
        """Updates Matplotlib plots, is called when stuff changed in Data"""
        if self.data is not None:
            self.ax.clear()
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
        self.ax.scatter(self.data[:, 0], self.data[:, 1], marker=".", s=2, c=self.selected_data)

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
            contained = path_list[last].contains_points(self.data)
            if self.DEBUG:
                print("[DEBUG] SelectGUI.on_mpl_press()")
                print("Nbr of found paths:", len(path_list))
                print(contained[0], np.shape(contained), '\n', np.shape(self.data))
            # TODO: Maybe use the contained array directly als self.selected_data
            self.selected_data = np.empty(np.shape(contained), dtype=bool)
            for i in range(len(contained)):
                if contained[i]:
                    self.selected_data[i] = True
                    # self.selected_data = np.hstack(self.selected_data, self.data[i])  # -> not useful because processed data is not fixed to fasta contigs
                else:
                    self.selected_data[i] = False
            self.update_plot(highlighted_cont=path_list[last])  # , col='green')
            self.analyze_widget.update_selected(self.selected_data)
            self.set_selected_cnt()

    def analyze_selected(self):
        """Takes Selected Datapoints and checks in MG Data for MG's and calculates coverage and contamination"""
        self.analyze_widget.update_selected(self.selected_data)
        if not self.analyze_widget.isVisible():
            self.analyze_widget.show()
        # TODO: Take Bool(0,1) Array self.selected_data and get the corresponding Fasta contigs, then apply found mg's and calculate stats
        # TODO change a bit

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        super().closeEvent(a0)
        self.cfg.on_close()
