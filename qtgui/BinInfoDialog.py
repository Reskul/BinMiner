from PyQt5.QtWidgets import *


class BinInfoDialog(QDialog):
    def __init__(self, parent: QMainWindow):
        super().__init__(parent)
        self.setWindowTitle("Selected Bin")
