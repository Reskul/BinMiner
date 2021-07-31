import sys
import platform
from qtgui import *
from PyQt5.QtWidgets import *

DEFAULT_WIDTH = 1280
DEFAULT_HEIGHT = 720

if __name__ == '__main__':
    app = QApplication(sys.argv)

    os = platform.system()

    size = app.primaryScreen().size()
    x = (size.width() / 2) - (DEFAULT_WIDTH / 2)
    y = (size.height() / 2) - (DEFAULT_HEIGHT / 2)
    mw = Window(x, y, DEFAULT_WIDTH, DEFAULT_HEIGHT, os)
    mw.show()

    sys.exit(app.exec())
