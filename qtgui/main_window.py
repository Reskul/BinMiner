from PyQt5.QtCore import Qt, pyqtSlot, QThreadPool

from .widgets import *
from lib import *


class ControllingWindow(QMainWindow):
    TAG = "ControllingWindow"
    STATUS_INPUT = 0
    STATUS_SELECT = 1
    STATUS = STATUS_INPUT

    def __init__(self, x: int, y: int, w: int, h: int, cfg: Configurator = None, parent=None,
                 debug=False, test=False):
        super().__init__(parent)
        # GENERAL Settings ----------
        self.DEBUG = debug
        self.TEST = test
        if cfg:
            self.cfg = cfg

        self.setWindowTitle("Sequence Mining Tool")
        self.setGeometry(x, y, w, h)

        self.contigs = None
        self.mgs = None

        self.select_widget = None
        # self.rating_dialog = BinInfoDialog()

        # Loading Screen ----------
        # Custom Widget by https://github.com/snowwlex/QtWaitingSpinner/blob/master/README.md
        # Py-version: https://github.com/z3ntu/QtWaitingSpinner
        self.loading_spinner = QtWaitingSpinner(self, centerOnParent=True, disableParentWhenSpinning=True)

        self.determine_widget()
        # self.create_menubar()

    def determine_widget(self):
        if self.STATUS == self.STATUS_INPUT:
            self.setCentralWidget(InputGUI(parent=self, cfg=self.cfg))
        elif self.STATUS == self.STATUS_SELECT:
            self.setCentralWidget(self.select_widget)

    def process_input(self, contig_path, coverage_path, kmere_path, perplexity, plotstate, fetchmg_respath=None,
                      prodigal_path=None, fetchmg_path=None):
        if fetchmg_respath is None and fetchmg_path is None:
            print(f"[ERROR] FetchMG Results or path to FetchMG Bin must be provided.")
        elif fetchmg_respath is not None:
            runnable = DataLoadingRunnable(self, contig_path, coverage_path, kmere_path, self.cfg.homepath, perplexity,
                                           plotstate, fetchmg_respath, debug=self.DEBUG)
            QThreadPool.globalInstance().start(runnable)
        elif fetchmg_path is not None:
            print(f"[ERROR] Not finished this part yet ;).")

        self.loading_spinner.start()

    @pyqtSlot(np.ndarray, np.ndarray, np.ndarray)
    def data_ready(self, contigs, mgs, datapoints):
        self.loading_spinner.stop()
        self.contigs = contigs
        self.mgs = mgs
        self.STATUS = self.STATUS_SELECT
        self.select_widget = SelectGUI(datapoints, contigs, mgs, parent=self, debug=self.DEBUG)
        self.determine_widget()

    @pyqtSlot(str)
    def data_failed(self, message):
        print("Loading Data Failed:", message)
