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
        self.data_path_le = QLineEdit()
        self.data_path_le.setReadOnly(True)
        self.read_in_layout = QHBoxLayout()
        self.read_in_layout.addWidget(file_lbl_1)
        self.read_in_layout.addWidget(self.data_path_le)
        self.read_in_layout.addWidget(filebrowser_btn_1)

        self.fasta_path_le = QLineEdit()
        self.fasta_path_le.setReadOnly(True)

        self.fasta_in_layout = QHBoxLayout()
        self.fasta_in_layout.addWidget(file_lbl_0)
        self.fasta_in_layout.addWidget(self.fasta_path_le)
        self.fasta_in_layout.addWidget(filebrowser_btn_0)

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
            header = reader.read_header_only(open(path, 'r'))
            contigs = [c.contig_pure for c in header]
            mg = MarkerGene(f.split('.')[0])
            mg.add_contigs(contigs)
            mgs.append(mg)

        # if self.DEBUG:
        #     print(f"[DEBUG] MG[0]:\n{mgs[0]}")

        # READ-IN Data from Fasta to get all existing Contigs
        reader = FastaReader(FastaReader.MYCC)
        header = reader.read_header_only(open(fasta_path, 'r'))
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
