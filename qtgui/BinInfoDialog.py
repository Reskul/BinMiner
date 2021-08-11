from PyQt5.QtWidgets import *


class BinInfoDialog(QDialog):
    def __init__(self, parent: QMainWindow, selected, debug=False):
        super().__init__(parent)
        self.setWindowTitle("Selected Bin")
        self.DEBUG = False
        self.data = selected

        self.update_gui()

    def update_selected(self, selected):
        self.data = selected
        self.update_gui()

    def update_gui(self):
        pass
