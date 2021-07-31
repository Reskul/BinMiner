from PyQt5 import QtCore
from PyQt5.QtWidgets import *


class PathsDialog(QDialog):
    OS = None
    PARENT = None
    DEFAULTPATH = None
    DEBUG = None

    def __init__(self, parent: QMainWindow, defaultpath, os: str = 'Linux', debug=False) -> None:
        super().__init__(parent)
        self.setModal(True)
        self.setWindowTitle("Set Paths")

        # TODO: Read maybe existing Paths from file and write to LE's

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

        # Control btns
        ok_btn = QPushButton("Ok")
        ok_btn.clicked.connect(self.ok_clicked)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.cancel_clicked)

        # Overall Layout ----------
        layout = QGridLayout(self)
        layout.addWidget(prodigal_lbl, 0, 0)
        layout.addWidget(self.prodigal_le, 0, 1)
        layout.addWidget(prodigal_btn, 0, 2)
        layout.addWidget(fetchmg_lbl, 1, 0)
        layout.addWidget(self.fetchMG_le, 1, 1)
        layout.addWidget(fetchmg_btn, 1, 2)
        layout.addWidget(data_lbl, 2, 0)
        layout.addWidget(self.data_le, 2, 1)
        layout.addWidget(data_btn, 2, 2)
        layout.addWidget(ok_btn, 3, 2)
        layout.addWidget(cancel_btn, 3, 0)

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

    def set_fetchmg(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select FetchMG Script', self.DEFAULTPATH,
                                              'Perl Scripts (*.pl)')
        if path:
            self.fetchMG_le.setText(path)
            self.PARENT.fetchMG_path = path

    def set_datadir(self):
        """select dir where all data is stored"""
        dir_path = None
        if self.OS == 'Linux':
            dir_path = QFileDialog.getExistingDirectory(self, 'Select Directory', '/home', QFileDialog.ShowDirsOnly)
        elif self.OS == 'Windows':
            dir_path = QFileDialog.getExistingDirectory(self, 'Select Directory', 'C:\\Users', QFileDialog.ShowDirsOnly)
        if dir_path:
            self.data_le.setText(dir_path)
            # self.PARENT.DATADIR = dir_path  # TODO: add this to main Window

    def ok_clicked(self):
        # TODO: Write paths to some file
        self.close()

    def cancel_clicked(self):
        self.close()
