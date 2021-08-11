import os
import subprocess
import numpy as np
from PyQt5 import QtGui
from PyQt5.QtCore import QMetaObject, Q_ARG, Qt, QRunnable, pyqtSlot, QThreadPool
from PyQt5.QtWidgets import *
from matplotlib import patches
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from lib import *
from sklearn.manifold import TSNE
from KDEpy import FFTKDE
from skimage.feature import peak_local_max
from .dialogs import *
from cfg import Configurator
import ntpath
from fasta import *


class Window(QMainWindow):
    DEFAULTPATH_WIN = "C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Data\\precomputed\\"  # TODO: remove kind of
    DEFAULTPATH_LINUX = "~/Dokumente/BA_Data/"
    DEFAULTPATH = None
    DEBUG = True
    OS = None

    DATADIR = None

    prodigal_path = None
    fetchMG_path = None

    TMP_fetchMG_results_path = None
    TMP_prodigal_results_path = None

    main_widget = None
    analyze_widget = None

    data_path = None
    data = None

    CONTIGS = None

    selected_data = None
    selected_nbr = 0
    grid_points = 100

    cont_data = None

    bw = 1
    contours = 25

    def __init__(self, x: int, y: int, w: int, h: int, os: str, cfg: Configurator = None, parent=None):
        super().__init__(parent)
        if cfg:
            self.cfg = cfg
        if os == 'Linux':
            self.DEFAULTPATH = self.DEFAULTPATH_LINUX
            if self.cfg:
                self.prodigal_path = self.cfg.read(self.cfg.PRODIGAL_KEY)
                self.fetchMG_path = self.cfg.read(self.cfg.FETCHMG_KEY)
                # self.DATADIR = self.cfg.read(self.cfg.DATA_KEY)
                self.DATADIR = self.cfg.homepath

        elif os == 'Windows':
            self.DEFAULTPATH = self.DEFAULTPATH_WIN
        self.OS = os

        self.is_processing_fasta = False

        # General Settings
        self.main_widget = QWidget()
        self.setWindowTitle("Sequence Mining Tool")
        self.setGeometry(x, y, w, h)
        self.setCentralWidget(self.main_widget)
        self.create_menubar()

        # Data Read-In ----------
        file_lbl_0 = QLabel("Fasta Datensatz wählen:")
        filebrowser_btn_0 = QPushButton("Datei wählen")
        filebrowser_btn_0.clicked.connect(self.select_fasta_file)
        file_lbl_1 = QLabel("Pre-computed Datensatz wählen:")
        filebrowser_btn_1 = QPushButton("Datei wählen")
        filebrowser_btn_1.clicked.connect(self.select_npy_file)

        self.fasta_path_le = QLineEdit()
        self.fasta_path_le.setReadOnly(True)
        self.data_path_le = QLineEdit()
        self.data_path_le.setReadOnly(True)

        self.fasta_in_layout = QHBoxLayout()
        self.fasta_in_layout.addWidget(file_lbl_0)
        self.fasta_in_layout.addWidget(self.fasta_path_le)
        self.fasta_in_layout.addWidget(filebrowser_btn_0)
        self.read_in_layout = QHBoxLayout()
        self.read_in_layout.addWidget(file_lbl_1)
        self.read_in_layout.addWidget(self.data_path_le)
        self.read_in_layout.addWidget(filebrowser_btn_1)

        # Data Visualization ----------
        fig = Figure()
        self.ax = fig.add_subplot(111)
        self.canvas = FigureCanvas(fig)
        # self.canvas.resize(300, 300)
        toolbar = NavigationToolbar(self.canvas, self.main_widget)

        diagram_layout = QVBoxLayout()
        diagram_layout.addWidget(self.canvas)
        diagram_layout.addWidget(toolbar)
        # Matplotlib Interaction ----------
        self.canvas.mpl_connect("button_press_event", self.on_mpl_press)

        # Loading Screen ----------
        # Custom Widget by https://github.com/snowwlex/QtWaitingSpinner/blob/master/README.md
        # Py-version: https://github.com/z3ntu/QtWaitingSpinner
        self.loading_spinner1 = QtWaitingSpinner(self, centerOnParent=False, disableParentWhenSpinning=True)
        self.loading_spinner2 = QtWaitingSpinner(self, centerOnParent=False, disableParentWhenSpinning=False)
        self.spacer_0 = QSpacerItem(self.loading_spinner1.sizeHint().width(), self.loading_spinner1.sizeHint().height())
        self.spacer_1 = QSpacerItem(self.loading_spinner2.sizeHint().width(), self.loading_spinner2.sizeHint().height())

        self.fasta_in_layout.addWidget(self.loading_spinner1)
        self.fasta_in_layout.addSpacerItem(self.spacer_0)

        self.read_in_layout.addWidget(self.loading_spinner2)
        self.read_in_layout.addSpacerItem(self.spacer_1)

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

        # Process Selected Layout ----------
        process_sel_btn = QPushButton("Ausgewählte Analysieren")
        process_sel_btn.clicked.connect(self.analyze_selected)
        self.sel_lbl = QLabel()
        self.set_selected_cnt()

        # Manage layout ----------
        layout = QGridLayout()
        layout.addLayout(self.fasta_in_layout, 0, 0, 1, 2)
        layout.addLayout(self.read_in_layout, 1, 0, 1, 2)
        layout.addLayout(diagram_layout, 2, 0)
        layout.addLayout(slider_layout, 2, 1)
        layout.addWidget(self.sel_lbl, 3, 0, alignment=Qt.AlignRight)
        layout.addWidget(process_sel_btn, 3, 1)

        self.main_widget.setLayout(layout)

    def set_selected_cnt(self):
        if self.selected_data is not None:
            count = sum(self.selected_data)
            self.selected_nbr = count
            self.sel_lbl.setText(f"Ausgewählt: {count}")
        else:
            self.sel_lbl.setText(f"Ausgewählt: --")

    def select_npy_file(self):
        # TODO Change to "Load k-Mer Data" or something like this
        path, _ = QFileDialog.getOpenFileName(self.main_widget, 'Open Numpy File', self.DEFAULTPATH,
                                              'Numpy Files (*.npy)')
        if self.DEBUG:
            print(path)
        if path:
            # print(_)
            self.data_path = path
            # This could take several seconds --> Loading Symbol is shown and UI deactivated
            self.data_path_le.setText(path)
            self.read_in_layout.removeItem(self.spacer_1)
            self.loading_spinner2.start()

            runnable = LoadingNpyRunnable(path, self, self.DEBUG)
            QThreadPool.globalInstance().start(runnable)

    def select_fasta_file(self):
        """Select Fasta File for fundamental Data input and start pipe"""
        # TODO: change to "Start Marker Gene calculation"
        if not self.is_processing_fasta:
            fasta_path, _ = QFileDialog.getOpenFileName(self.main_widget, 'Open Fasta File', self.DEFAULTPATH,
                                                        'Fasta Files (*.fasta)')
            self.fasta_path_le.setText(fasta_path)

            runnable = FastaLoadingRunnable(fasta_path, self, self.prodigal_path, self.fetchMG_path, self.DATADIR,
                                            self.DEBUG)
            QThreadPool.globalInstance().start(runnable)

    def create_menubar(self):
        """Creates Actions and the Menubar to show them in"""
        set_paths_act = QAction("Set Paths", self)
        set_paths_act.triggered.connect(self.start_paths_dialog)

        # use_existing_mgs = QAction("Use Existing Files", self)
        # use_existing_mgs.triggered.connect(self.process_markergenes)

        menubar = self.menuBar()
        settings_menu = QMenu("&Settings", self)
        settings_menu.addAction(set_paths_act)
        # settings_menu.addAction(use_existing_mgs)

        menubar.addMenu(settings_menu)
        help_menu = menubar.addMenu("Help")
        # TODO fill Help Menu

    def start_paths_dialog(self):
        dialog = PathsDialog(self, self.DEFAULTPATH, self.cfg, debug=self.DEBUG)
        dialog.show()

    @pyqtSlot(np.ndarray)
    def set_data(self, data):
        if self.DEBUG:
            print("[DEBUG] Window.set_data()")
        self.data = data
        self.selected_data = np.zeros(len(data))
        self.loading_spinner2.stop()
        self.read_in_layout.addSpacerItem(self.spacer_1)
        self.update_plot()

    @pyqtSlot(np.ndarray)
    def protdata_ready(self, contigs):
        if self.DEBUG:
            print("[DEBUG] Window.protdata_ready()")
        self.loading_spinner1.stop()
        self.fasta_in_layout.addSpacerItem(self.spacer_0)
        # TODO: show grüner haken lol
        self.CONTIGS = contigs
        self.is_processing_fasta = False

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        super().closeEvent(a0)
        self.cfg.on_close()

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
        scatter_data = self.ax.scatter(self.data[:, 0], self.data[:, 1], marker=".", s=2, c=self.selected_data)

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
                print("Nbr of found paths:", len(path_list))
                print(contained[0], np.shape(contained), '\n', np.shape(self.data))
            # TODO: Maybe use the contained array directly als self.selected_data
            self.selected_data = np.empty(np.shape(contained))
            for i in range(len(contained)):
                if contained[i]:
                    self.selected_data[i] = 1
                    # self.selected_data = np.hstack(self.selected_data, self.data[i])  # -> not useful because processed data is not fixed to fasta contigs
                else:
                    self.selected_data[i] = 0
            self.update_plot(highlighted_cont=path_list[last])  # , col='green')
            self.set_selected_cnt()

    def process_markergenes(self, tmp_mg_res: str, tmp_fasta_path: str, tmp_fasta_translation_path: str):
        self.is_processing_fasta = True
        self.fasta_in_layout.removeItem(self.spacer_0)
        self.loading_spinner1.start()
        translate = False
        if tmp_fasta_translation_path:
            translate = True
        runnable = FastaLoadingRunnable(tmp_fasta_path, self, self.prodigal_path, tmp_mg_res, self.DATADIR,
                                        debug=self.DEBUG, only_analyze=True, contig_translation=translate,
                                        translation_file=tmp_fasta_translation_path)
        QThreadPool.globalInstance().start(runnable)

    def analyze_selected(self):
        """Takes Selected Datapoints and checks in MG Data for MG's and calculates coverage and contamination"""
        self.analyze_widget = BinInfoDialog(self, self.selected_data, self.DEBUG)
        self.analyze_widget.show()
        # TODO: Take Bool(0,1) Array self.selected_data and get the corresponding Fasta contigs, then apply found mg's and calculate stats
        pass


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
