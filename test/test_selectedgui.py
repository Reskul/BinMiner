from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow
import numpy as np
from qtgui import *
import sys
import os
from sklearn.manifold import TSNE


class DummyMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        path = "~/Dokuments/BA_Arbeit/"
        path = os.path.expanduser(path)

        path, _ = QFileDialog.getOpenFileName(self, 'Select Precomputed K-mere Data', path,
                                              'Fasta Files (*.fasta')
        if path:
            data_raw = np.load(path)
            contig_lengths = np.sum(data_raw, 1)
            x_mat = (data_raw.T / contig_lengths).T  # Norm data ??? ASK!
            data = TSNE(n_components=2).fit_transform(x_mat)  # T-SNE Data Dim reduction
            wdgt = SelectGUI(data)
            self.setCentralWidget(wdgt)
        else:
            print("Big fail Junge!")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    mw = DummyMainWindow()
    mw.show()
    sys.exit(app.exec())
