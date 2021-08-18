import numpy as np
import os
from PyQt5 import QtGui
from PyQt5.QtCore import Qt, pyqtSlot, QThreadPool
from matplotlib import patches
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from KDEpy import FFTKDE
from skimage.feature import peak_local_max

from lib import *
from .dialogs import *
from cfg import Configurator
from worker import *


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

    data_path = None
    data = None
    access = None

    CONTIGS = None
    MGS = None

    selected_data = None
    selected_nbr = 0
    grid_points = 100

    cont_data = None

    bw = 1
    contours = 25

    def __init__(self, x: int, y: int, w: int, h: int, operating_system: str, cfg: Configurator = None, parent=None):
        super().__init__(parent)
        if cfg:
            self.cfg = cfg
        if operating_system == 'Linux':
            self.DEFAULTPATH = os.path.expanduser(self.DEFAULTPATH_LINUX)
            if self.cfg:
                self.prodigal_path = self.cfg.read(self.cfg.PRODIGAL_KEY)
                self.fetchMG_path = self.cfg.read(self.cfg.FETCHMG_KEY)
                # self.DATADIR = self.cfg.read(self.cfg.DATA_KEY)
                self.DATADIR = self.cfg.homepath

        elif operating_system == 'Windows':
            self.DEFAULTPATH = self.DEFAULTPATH_WIN
        self.OS = operating_system

        self.analyze_widget = BinInfoDialog(self, debug=self.DEBUG)

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

    @pyqtSlot(np.ndarray, np.ndarray)
    def set_data(self, data, access):
        if self.DEBUG:
            print("[DEBUG] Window.set_data()")
        self.data = data
        self.access = access
        self.selected_data = np.zeros(len(data))
        self.loading_spinner2.stop()
        self.read_in_layout.addSpacerItem(self.spacer_1)
        self.update_plot()

    @pyqtSlot(np.ndarray, np.ndarray)
    def protdata_ready(self, contigs, mgs):
        if self.DEBUG:
            print("[DEBUG] Window.protdata_ready()")
        self.loading_spinner1.stop()
        self.fasta_in_layout.addSpacerItem(self.spacer_0)
        # TODO: show grüner haken lol
        self.CONTIGS = contigs
        self.MGS = mgs
        self.analyze_widget.set_contigs_and_markergenes(contigs, mgs)
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
        self.analyze_widget.update_selected(self.selected_data)
        if not self.analyze_widget.isVisible():
            self.analyze_widget.show()
        # TODO: Take Bool(0,1) Array self.selected_data and get the corresponding Fasta contigs, then apply found mg's and calculate stats
        pass


class ControllingWindow(QMainWindow):
    MAINWIDGET = None

    def __init__(self):
        pass
