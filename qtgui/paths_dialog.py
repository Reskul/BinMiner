from PyQt5.QtWidgets import *
from cfg import *


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
        self.setFixedWidth(720)
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
        mg_btn.clicked.connect(self.find_fetchmg_results)
        self.mg_le = QLineEdit()
        self.mg_le.setReadOnly(True)

        fasta_lbl = QLabel("Original Daten:")
        fasta_btn = QPushButton("Wählen")
        fasta_btn.clicked.connect(self.find_input_data)
        self.fasta_le = QLineEdit()
        self.fasta_le.setReadOnly(True)

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
        layout.addWidget(start_btn, 6, 2)

        layout.addWidget(ok_btn, 7, 2)
        layout.addWidget(cancel_btn, 7, 0)

        # Write Text to LE's if path is already set ----------
        # self.prodigal_le.setText(self.cfg.read(self.cfg.PRODIGAL_KEY))
        # self.fetchMG_le.setText(self.cfg.read(self.cfg.FETCHMG_KEY))
        # self.data_le.setText(self.cfg.read(self.cfg.DATA_KEY))

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

    def find_fetchmg_results(self):
        path = QFileDialog.getExistingDirectory(self, 'Select Directory', self.cfg.homepath, QFileDialog.ShowDirsOnly)
        if path:
            self.mg_le.setText(path)
            self.PARENT.TMP_fetchMG_results_path = path

    def find_input_data(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select Fasta File', self.cfg.homepath, "Fasta Files (*.fasta)")
        if path:
            self.fasta_le.setText(path)
            self.PARENT.TMP_prodigal_results_path = path

    def start_result_processing(self):
        self.PARENT.process_markergenes(self.mg_le.text(), self.fasta_le.text())

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
