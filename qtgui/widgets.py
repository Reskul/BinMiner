from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from matplotlib import patches
from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.backend_bases import _Mode
from matplotlib.figure import Figure
import numpy as np

from KDEpy import FFTKDE
from KDEpy.BaseKDE import BaseKDE
from KDEpy.bw_selection import improved_sheather_jones, silvermans_rule
from skimage.feature import peak_local_max


from lib import *
from .dialogs import BinInfoDialog, NameSelectedDialog


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


class MyNavigationToolbar(NavigationToolbar2QT):
    def _init_toolbar(self):
        pass

    def home(self, *args):
        super().home(args)
        self.canvas.parent().x_lim = self.canvas.parent().home_x_lim
        self.canvas.parent().y_lim = self.canvas.parent().home_y_lim


class InputGUI(QWidget):
    PLOTSTATE_KMERE = 0
    PLOTSTATE_COV = 1
    PLOTSTATE_COMB = 2

    def __init__(self, parent, cfg: Configurator = None, debug=False):
        super().__init__(parent)
        # INSTANCE VARIABLES
        self.DEBUG = debug
        self.TEST = parent.TEST
        self.PLOTSTATE = 0
        self.config = cfg
        self.last_path = self.config.homepath
        self.contig_sequences_path = None
        self.contig_coverage_path = None
        self.kmere_dataset_path = None
        self.test_path = None
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
        f_layout.addRow(QLabel("Contig sequence File (.fasta)"), self.contig_dataset_le)

        # Contig Coverage File
        self.contig_coverage_le = QFileInputLine()
        self.contig_coverage_le.clicked.connect(self.cc_clicked)
        f_layout.addRow((QLabel("Contig coverage file (.txt)")), self.contig_coverage_le)

        if self.DEBUG:
            # K-Mere Data ----------
            self.kmere_dataset_le = QFileInputLine()
            self.kmere_dataset_le.clicked.connect(self.km_clicked)
            f_layout.addRow(QLabel("K-mer data file (.npy)"), self.kmere_dataset_le)

        # Contig Test Data
        if self.TEST:
            self.test_le = QFileInputLine()
            self.test_le.clicked.connect(self.test_clicked)
            f_layout.addRow(QLabel("Contig labels file (.txt)"), self.test_le)

        # TSNE Perplexity Parameter
        perplex_lbl = QLabel("T-SNE perplexity")
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
        self.source_radbtn = QRadioButton("Calculate MarkerGenes(linux only)")
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

        # SELECT METHOD
        method_gbox = QGroupBox("Method")
        kmere_plot_radbtn = QRadioButton("Plot K-mer data")
        kmere_plot_radbtn.setChecked(True)
        cov_plot_radbtn = QRadioButton("Plot coverage data")
        combined_plot_radbtn = QRadioButton("Plot combined data")

        kmere_plot_radbtn.clicked.connect(self.radio_kmere_clicked)
        cov_plot_radbtn.clicked.connect(self.radio_cov_clicked)
        combined_plot_radbtn.clicked.connect(self.radio_comb_clicked)

        method_vbox = QVBoxLayout()
        method_vbox.addWidget(kmere_plot_radbtn)
        method_vbox.addWidget(cov_plot_radbtn)
        method_vbox.addWidget(combined_plot_radbtn)

        method_gbox.setLayout(method_vbox)

        # "Next" Button
        next_btn = QPushButton("Next")
        next_btn.clicked.connect(self.next_clicked)

        # Same Procedure as last year Miss Sophie?
        self.last_ckbox = QCheckBox("Use last values")
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
        general_layout.addWidget(method_gbox)
        general_layout.addLayout(btn_layout)
        general_layout.addStretch(1)
        self.setLayout(general_layout)

    def cd_clicked(self):
        # path = QFileDialog.getExistingDirectory(self, 'Select Contig Fasta-Data File', self.config.homepath,
        #                                        QFileDialog.ShowDirsOnly)
        path, _ = QFileDialog.getOpenFileName(self, 'Select Contig fasta-data file', self.last_path,
                                              'Fasta Files (*.fasta *.faa *.fa *.fna)')
        if path:
            self.contig_dataset_le.setText(path)
            self.contig_sequences_path = path
            self.last_path = path

    def cc_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Contig coverage file', self.last_path,
                                              'Depth/Text Files (*.depth *.depth.txt)')
        if path:
            self.contig_coverage_le.setText(path)
            self.contig_coverage_path = path
            self.last_path = path

    def km_clicked(self):
        # path = QFileDialog.getExistingDirectory(self, 'Select K-mere Data Directory', self.config.homepath,
        #                                         QFileDialog.ShowDirsOnly)
        path, _ = QFileDialog.getOpenFileName(self, 'Select precomputed K-mer data', self.last_path,
                                              'Numpy Files (*.npy)')
        if path:
            self.kmere_dataset_le.setText(path)
            self.kmere_dataset_path = path
            self.last_path = path

    def test_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Contig Test Data', self.last_path,
                                              'Text Files (*.txt)')
        if path:
            self.test_le.setText(path)
            self.test_path = path
            self.last_path = path

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
        path = QFileDialog.getExistingDirectory(self, 'Select FetchMG result directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
        if path:
            self.fetchmg_res_le.setText(path)
            self.fetchmg_result_path = path

    def pp_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Prodigal binary', self.DEFAULTPATH,
                                              'Binaries (*.linux)')
        if path:
            self.prodigal_binpath = path
            self.prodigal_path_le.setText(path)

    def fmgp_clicked(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select FetchMG script', self.DEFAULTPATH,
                                              'Perl Scripts (*.pl)')
        if path:
            self.fetchmg_binpath = path
            self.fetchmg_path_le.setText(path)

    def radio_kmere_clicked(self):
        self.PLOTSTATE = self.PLOTSTATE_KMERE

    def radio_cov_clicked(self):
        self.PLOTSTATE = self.PLOTSTATE_COV

    def radio_comb_clicked(self):
        self.PLOTSTATE = self.PLOTSTATE_COMB

    def next_clicked(self):
        if self.contig_sequences_path is not None:
            if self.results_radbtn.isChecked() and self.fetchmg_result_path is not None:
                if self.config:
                    self.config.write(self.config.CONTIGSEQ_KEY, self.contig_sequences_path)
                    self.config.write(self.config.CONTIGCOV_KEY, self.contig_coverage_path)
                    self.config.write(self.config.KMERE_KEY, self.kmere_dataset_path)
                    self.config.write(self.config.FETCHMGRES_KEY, self.fetchmg_result_path)
                self.parent().process_input(self.contig_sequences_path, self.contig_coverage_path,
                                            self.perplex_spbox.value(), self.PLOTSTATE, kmere_path=self.kmere_dataset_path,
                                            fetchmg_respath=self.fetchmg_result_path, testdata_path=self.test_path)

            elif self.source_radbtn.isChecked() and self.fetchmg_binpath is not None and self.prodigal_binpath is not None:
                if self.cfg:
                    self.config.write(self.config.CONTIGSEQ_KEY, self.contig_sequences_path)
                    self.config.write(self.config.CONTIGCOV_KEY, self.contig_coverage_path)
                    self.config.write(self.config.KMERE_KEY, self.kmere_dataset_path)
                    self.config.write(self.config.FETCHMG_KEY, self.fetchmg_binpath)
                    self.config.write(self.config.PRODIGAL_KEY, self.prodigal_binpath)
                self.parent().process_input(self.contig_sequences_path, self.contig_coverage_path,
                                            self.perplex_spbox.value(), self.PLOTSTATE, kmere_path=self.kmere_dataset_path,
                                            prodigal_path=self.prodigal_binpath, fetchmg_path=self.fetchmg_binpath, testdata_path=self.test_path)
            else:
                print(f"[ERROR] Something went terribly wrong with radio buttons.")

    def lastckbox_clicked(self):
        if self.last_ckbox.isChecked():
            # Fill LineEdits with Values
            if not self.config.is_new:
                if self.DEBUG:
                    print(f"[DEBUG] InputGUI.lastckbox_clicked(): Config wurde NICHT neu erstellt.")
                contig_seq = self.config.read(Configurator.CONTIGSEQ_KEY)
                contig_cov = self.config.read(Configurator.CONTIGCOV_KEY)
                # kmere_set = self.config.read(Configurator.KMERE_KEY)
                fetchmg_res = self.config.read(Configurator.FETCHMGRES_KEY)
                fetchmg_bin = self.config.read(Configurator.FETCHMG_KEY)
                prodigal_bin = self.config.read(Configurator.PRODIGAL_KEY)

                self.contig_sequences_path = contig_seq
                self.contig_coverage_path = contig_cov
                # self.kmere_dataset_path = kmere_set
                self.fetchmg_result_path = fetchmg_res
                self.fetchmg_binpath = fetchmg_bin
                self.prodigal_binpath = prodigal_bin

                self.contig_dataset_le.setText(contig_seq)
                self.contig_coverage_le.setText(contig_cov)
                # self.kmere_dataset_le.setText(kmere_set)
                self.fetchmg_res_le.setText(fetchmg_res)
                self.fetchmg_path_le.setText(fetchmg_bin)
                self.prodigal_path_le.setText(prodigal_bin)
        else:
            # clear all LineEdits
            self.contig_dataset_le.setText("")
            self.contig_coverage_le.setText("")
            # self.kmere_dataset_le.setText("")
            self.fetchmg_res_le.setText("")
            self.fetchmg_path_le.setText("")
            self.prodigal_path_le.setText("")

            self.contig_sequences_path = None
            self.contig_coverage_path = None
            # self.kmere_dataset_path = None
            self.fetchmg_result_path = None
            self.fetchmg_binpath = None
            self.prodigal_binpath = None


class SelectGUI(QWidget):
    def __init__(self, kmere_data, contigs, mgs, parent=None, debug=False, test=False):
        super().__init__(parent)
        # VARIABLE Initialization
        self.DEBUG = debug
        self.TEST = test
        self.data = kmere_data
        self.contigs = contigs
        self.covdim = len(contigs[0].coverage)
        self.mgs = mgs
        self.mg_dict = {}
        for i_idx in range(len(mgs)):
            self.mg_dict[mgs[i_idx].MG_name] = i_idx
        self.cont_data = None
        self.n_selected = 0
        self.completeness = None
        self.contamination = None
        self.contours = 25
        self.selected_vec = None
        self.grid_points = 100

        self.contour_visible = False
        self.patch = None

        self.prototype_name = None

        self.x_lim = None
        self.y_lim = None

        self.analyze_widget = BinInfoDialog(self.parent(), debug=self.DEBUG)

        # Calculating Color and size values for depiction, Color-> Coverage & Size-> Contig length
        self.colormap = np.empty(len(contigs))
        covs = np.empty(len(contigs))
        self.sizemap = np.empty(len(contigs))
        i = 0
        for contig in self.contigs:
            self.sizemap[i] = np.log2(len(contig.sequence))  # TODO: Work out a good scale
            covs[i] = contig.coverage_1d
            i += 1

        if self.DEBUG:
            print(f"[DEBUG] SelectGUI.__init__(): Log-Length's of sequences:")
            print(self.sizemap)

        covmin = np.min(covs)
        covmax = np.max(covs)
        covdiff = covmax - covmin
        for i in range(len(contigs)):
            self.colormap[i] = (covs[i] - covmin) / covdiff

        if self.DEBUG:
            print(self.colormap)

        # Bandwith calculation
        prep_data_a = BaseKDE._process_sequence(self.data[:, 0])
        prep_data_b = BaseKDE._process_sequence(self.data[:, 1])

        self.bw_lower_bound = min(improved_sheather_jones(prep_data_a), improved_sheather_jones(prep_data_b))
        self.bw_upper_bound = max(silvermans_rule(prep_data_a), silvermans_rule(prep_data_b))
        bw_diff = self.bw_upper_bound - self.bw_lower_bound

        if self.DEBUG:
            print(f"[DEBUG] SelectGUI.__init__(): LowerBW:{self.bw_lower_bound} | UpperBW:{self.bw_upper_bound} | Diff:{bw_diff}")

        self.bw = self.bw_lower_bound + bw_diff / 2
        n_bw_slider_steps = int(np.round(bw_diff / 0.1))

        # BACK BUTTON
        back_btn = QPushButton("<- Back")
        back_btn.clicked.connect(self.back_clicked)
        back_btn.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # Data Visualization ----------
        fig = Figure()
        self.ax: Axes = fig.add_subplot(111)
        self.canvas = FigureCanvas(fig)
        # self.canvas.resize(300, 300)
        self.toolbar = MyNavigationToolbar(self.canvas, self)

        diagram_layout = QGridLayout()
        diagram_layout.addWidget(back_btn, 0, 0, 1, 1)
        diagram_layout.addWidget(self.canvas, 1, 0, 1, 3)
        diagram_layout.addWidget(self.toolbar, 2, 0, 1, 1)

        # Matplotlib Interaction ----------
        self.canvas.mpl_connect("button_press_event", self.on_mpl_press)
        self.canvas.mpl_connect("button_release_event", self.on_mpl_release)

        # Info Bar Layout ----------
        process_sel_btn = QPushButton("Check selected")
        process_sel_btn.clicked.connect(self.show_diagram_dialog)
        self.sel_lbl = QLabel("Selected: --")
        self.completeness_lbl = QLabel("Completeness: --")
        self.contamination_lbl = QLabel("Contamination: --")

        info_bar_layout = QHBoxLayout()
        info_bar_layout.addWidget(self.completeness_lbl)
        info_bar_layout.addWidget(self.contamination_lbl)
        info_bar_layout.addWidget(self.sel_lbl)
        info_bar_layout.addWidget(process_sel_btn)

        diagram_layout.addLayout(info_bar_layout, 2, 2, 1, 1)
        self.update_info_bar()

        # Save Selection to fasta btn (in Slider Box) ---------
        save_selected_btn = QPushButton("save selected")
        save_selected_btn.clicked.connect(self.save_pressed)
        # Slider Section ----------
        bw_slider_lbl = QLabel("Bandwidth")
        cont_slider_lbl = QLabel("Contours")
        self.bw_nbr_lbl = QLabel("{:10.2f}".format(self.bw))
        self.cont_nbr_lbl = QLabel(str(self.contours))
        self.bw_slider = QSlider(Qt.Vertical)
        self.cont_slider = QSlider(Qt.Vertical)

        self.bw_slider.setRange(1, n_bw_slider_steps)
        self.bw_slider.setValue(int(n_bw_slider_steps / 2))
        self.cont_slider.setRange(1, 50)
        self.cont_slider.setValue(25)

        # self.bw_slider.sliderReleased.connect(self.on_bw_slider_change)
        self.bw_slider.sliderMoved.connect(self.on_bw_slider_change)
        # self.cont_slider.sliderReleased.connect(self.on_cont_slider_change)
        self.cont_slider.sliderMoved.connect(self.on_cont_slider_change)

        slider_layout = QGridLayout()  # TODO Add Button underneath to print selected to fasta file, maybe also somehow allow naming the "found" organism
        slider_layout.addWidget(self.bw_nbr_lbl, 0, 0, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(self.bw_slider, 1, 0, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(bw_slider_lbl, 2, 0)
        slider_layout.addWidget(self.cont_nbr_lbl, 0, 1, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(self.cont_slider, 1, 1, alignment=Qt.AlignHCenter)
        slider_layout.addWidget(cont_slider_lbl, 2, 1)
        slider_layout.addWidget(save_selected_btn, 3, 0, 1, 2)

        # GENERAL LAYOUT
        layout = QHBoxLayout()
        layout.addLayout(diagram_layout)
        layout.addLayout(slider_layout)

        self.update_plot()
        self.home_x_lim = self.ax.get_xlim()
        self.home_y_lim = self.ax.get_ylim()
        self.setLayout(layout)

    def update_info_bar(self):
        if self.selected_vec is not None:
            count = sum(self.selected_vec)
            self.n_selected = count
            self.sel_lbl.setText(f"Selected: {count}")
            self.completeness_lbl.setText(f"Completeness: {self.completeness}")
            self.contamination_lbl.setText(f"Contamination: {self.contamination}")
        else:
            self.sel_lbl.setText(f"Selected: --")
            self.completeness_lbl.setText("Completeness: --")
            self.contamination_lbl.setText("Contamination: --")

    def update_plot(self, selection: Path = None, col=None):
        """Updates Matplotlib plots, is called when stuff changed in Data"""
        if self.data is not None:
            self.ax.clear()
            if self.x_lim is not None and self.y_lim is not None:
                self.ax.set_xlim(left=self.x_lim[0], right=self.x_lim[1])
                self.ax.set_ylim(bottom=self.y_lim[0], top=self.y_lim[1])
                if self.DEBUG:
                    print(f"[DEBUG] SelectGUI.update_plot(): Setting Axes View Limits")
            self.ax.tick_params(axis='x', labelsize='14')
            self.ax.tick_params(axis='y', labelsize='14')

            # x, y, z = self.update_kde()
            # if highlighted_cont is not None:
            #     self.ax.patches.pop()
            #     patch = patches.PathPatch(highlighted_cont, facecolor=col, lw=1, edgecolor='black', fill=False)
            #     self.ax.add_patch(patch)
            # TODO find a way to draw a outline without convex hull, because convex hull could show the user wrong points that are not really selected!

            # update kernel density estimation contours
            x, y, z = self.update_kde()
            if self.DEBUG:
                print("[DEBUG] SelectGUI.update_plot(): KDE X:", x.shape, " Y:", y.shape, " Z:", z.shape)
            contour_coords = np.meshgrid(x, y)
            if selection is not None:
                if self.DEBUG:
                    print("[DEBUG]\tSelectGUI.update_plot(): Selection Path.codes:", selection.codes, "\n\t\tPath.vertices.shape:",
                          selection.vertices.shape)
                if self.contour_visible:
                    self.patch.remove()
                    self.contour_visible = False

                closepoly_idx = np.argwhere(selection.codes == 79)
                if self.DEBUG:
                    print("[DEBUG]\tSelectGUI.update_plot(): Closepoly Idx:", closepoly_idx)

                outline_vertices = selection.vertices[:closepoly_idx[0][0]+1]

                outline_path = Path(outline_vertices, codes=selection.codes[:closepoly_idx[0][0]+1], closed=True)
                self.patch = patches.PathPatch(outline_path, facecolor=col, lw=1, edgecolor='black', fill=False)
                self.ax.add_patch(self.patch)
                self.contour_visible = True

            # update the rest
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
        self.ax.scatter(self.data[:, 0], self.data[:, 1], marker=".", s=self.sizemap, c=self.colormap)

    def update_peaks(self, x, y, z):
        """Only updates the peak points"""
        peakdata = peak_local_max(z, threshold_rel=0.02)
        self.ax.scatter(x[peakdata[:, 1]], y[peakdata[:, 0]], marker=10, c="orange")

    def on_bw_slider_change(self):
        step = self.bw_slider.value()
        self.bw = self.bw_lower_bound + step * 0.1
        self.bw_nbr_lbl.setText("{:10.2f}".format(self.bw))
        self.update_plot()

    def on_cont_slider_change(self):
        self.contours = self.cont_slider.value()
        self.cont_nbr_lbl.setText(str(self.contours))
        self.update_plot()

    def on_mpl_press(self, e):
        """Matplotlib 'press' Event, calculates path and points selected by event"""
        if self.toolbar.mode == _Mode.NONE:
            if self.DEBUG:
                print(f"[DEBUG] MainWindow.on_mpl_press(): X:{e.xdata} Y:{e.ydata}")

            path_list = []
            if self.cont_data is not None:
                for path_collection in self.cont_data.collections:
                    check, _ = path_collection.contains(e)
                    if check:
                        # path.contains_point()->bool oder path.contains_points()-> bool array
                        paths = path_collection.get_paths()
                        for p in paths:
                            if p.contains_point((e.xdata, e.ydata)):
                                path_list.append(p)

                last = len(path_list) - 1
                self.selected_vec = np.array(path_list[last].contains_points(self.data), dtype=bool)
                if self.DEBUG:
                    print("[DEBUG] SelectGUI.on_mpl_press()")
                    print("Nbr of found paths:", len(path_list))
                    print(self.selected_vec[0], np.shape(self.selected_vec), '\n', np.shape(self.data))

                self.update_plot(selection=path_list[last])  # , col='green')
                sel_coverage, sel_kmer_counts = self.calc_values()
                self.analyze_widget.update_data(sel_kmer_counts, sel_coverage, self.covdim)
                self.update_info_bar()

    def on_mpl_release(self, e):
        if self.toolbar.mode == _Mode.ZOOM or self.toolbar.mode == _Mode.PAN:
            if self.DEBUG:
                print(f"[DEBUG] SelectGUI.on_mpl_release(): Getting Axes View Limits.")
            self.x_lim = self.ax.get_xlim()
            self.y_lim = self.ax.get_ylim()

    def calc_values(self):
        sel_contigs = self.contigs[self.selected_vec]
        sel_coverages = []
        sel_kmer_counts = []
        n_mgs = len(self.mgs)
        contained_mgs = np.zeros(n_mgs, dtype=int)  # counting which markergenes exist in selection
        for c in sel_contigs:
            c_mgs = c.mgs
            for mg in c_mgs:
                contained_mgs[self.mg_dict[mg]] += 1
            sel_coverages.append(c.coverage_1d)
            sel_kmer_counts.append(c.kmere_counts)

        val_greater_zero = [val > 0 for val in contained_mgs]
        self.completeness = sum(val_greater_zero) / n_mgs
        if self.DEBUG:
            print(f"[DEBUG] SelectGUI.calc_values()\n"
                  f"\tCounted MG's:{contained_mgs}\n"
                  f"\tValues greater than 0:{val_greater_zero}\n"
                  f"\tCompleteness:{self.completeness}")
        max_cnt = max(contained_mgs)
        i_max = 2
        contam = 0
        while i_max <= max_cnt:
            existing = sum([val == i_max for val in contained_mgs])
            contam += existing * (i_max - 1)
            i_max += 1
        self.contamination = contam / n_mgs

        if self.DEBUG:
            print(f"[DEBUG] SelectGUI.calc_values(): Contamination:{self.contamination}")

        return sel_coverages, sel_kmer_counts

    def save_pressed(self):
        dialog = NameSelectedDialog(self, self.DEBUG)
        dialog.show()

    def save_to_file(self, name='prototype', covdim=1):
        sel_contigs = self.contigs[self.selected_vec]
        n = len(sel_contigs)
        filepath, _ = QFileDialog.getSaveFileName(self, "Save selection to fasta file", f'{name}_n{n}.fasta')
        if filepath:
            if self.DEBUG:
                print(f"[VERBOSE] SelectGUI.save_selected_to_file(): Assembling & creating fasta file of selected contigs.")

            collection = []
            for contig in sel_contigs:
                if not self.TEST:
                    contig.organism = name
                collection.append(contig.to_fasta(covdim=covdim))
            content = "\n".join(collection).strip()
            file = open(filepath, 'w+')
            file.write(content)
            file.flush()
            file.close()

    def show_diagram_dialog(self):
        """Takes Selected Datapoints and checks in MG Data for MG's and calculates coverage and contamination"""
        if self.n_selected > 0:
            if not self.analyze_widget.isVisible():
                self.analyze_widget.show()

    def back_clicked(self):
        self.analyze_widget.close()
        self.parent().STATUS = self.parent().STATUS_INPUT
        self.parent().determine_widget()

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        super().closeEvent(a0)
        self.cfg.on_close()
