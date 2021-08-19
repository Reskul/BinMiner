from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from cfg import *


class QFileInputLine(QLineEdit):
    clicked = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)

    def mousePressEvent(self, a0: QtGui.QMouseEvent) -> None:
        self.clicked.emit()


class InputGUI(QWidget):
    def __init__(self, parent=None, cfg: Configurator = None):
        super().__init__(parent)
        # INSTANCE VARIABLES
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
        path = QFileDialog.getExistingDirectory(self, 'Select K-mere Data Directory', self.config.homepath,
                                                QFileDialog.ShowDirsOnly)
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
        self.parent().data_ready(self.contig_dataset_path, self.kmere_dataset_path)


class SelectGUI(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
