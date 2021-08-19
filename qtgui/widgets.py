from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from cfg import *
from .main_window import ControllingWindow


class InputGUI(QWidget):
    def __init__(self, parent: ControllingWindow = None, cfg: Configurator = None):
        super().__init__(parent)
        # INSTANCE VARIABLES
        self.config = cfg
        self.contig_dataset_path = None
        self.kmere_dataset_path = None

        # INIT GUI ---------- ---------
        f_layout = QFormLayout()
        form_gbox = QGroupBox()
        general_layout = QHBoxLayout()

        # Contig DNA Sequence ----------
        self.contig_dataset_le = QLineEdit()
        self.contig_dataset_le.setReadOnly(True)
        self.contig_dataset_le.clicked.connect(self.cd_clicked)
        f_layout.addRow(QLabel("Contig Dataset"), self.contig_dataset_le)

        # K-Mere Data ----------
        self.kmere_dataset_le = QLineEdit()
        self.kmere_dataset_le.setReadOnly(True)
        self.kmere_dataset_le.clicked.connect(self.km_clicked)
        f_layout.addRow(QLabel("K-Mere Dataset"), self.kmere_dataset_le)

        # "Next" Button
        next_btn = QPushButton("Next")
        next_btn.clicked.connect(self.next_clicked)

        # FINALIZE GUI
        form_gbox.setLayout(f_layout)
        general_layout.addWidget(form_gbox)
        general_layout.addWidget(next_btn, alignment=Qt.AlignRight)

    def cd_clicked(self):
        # path, _ = QFileDialog.getOpenFileName(self, 'Open Contig File', self.DEFAULTPATH,'Numpy Files (*.fasta)')
        path = QFileDialog.getExistingDirectory(self, 'Select Contig Data Directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
        self.contig_dataset_le.setText(path)
        self.contig_dataset_path = path

    def km_clicked(self):
        path = QFileDialog.getExistingDirectory(self, 'Select K-mere Data Directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
        self.kmere_dataset_le.setText(path)
        self.kmere_dataset_path = path

    def next_clicked(self):
        self.parent().data_ready(self.contig_dataset_path, self.kmere_dataset_path)


class SelectGUI(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
