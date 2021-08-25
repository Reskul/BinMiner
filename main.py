import sys
import os
import platform
from pathlib import Path
from lib import *
from qtgui import *
from PyQt5.QtWidgets import *

DEFAULT_WIDTH = 1280
DEFAULT_HEIGHT = 720

if __name__ == '__main__':
    app = QApplication(sys.argv)

    operating_system = platform.system()

    working_path = os.path.expanduser("~/Documents/MGB")
    if operating_system == 'Linux':
        if not os.path.exists(working_path):
            os.makedirs(working_path)
        else:
            print(f"{working_path} already exists.")

    cfg = Configurator(working_path)

    size = app.primaryScreen().size()
    x = (size.width() / 2) - (DEFAULT_WIDTH / 2)
    y = (size.height() / 2) - (DEFAULT_HEIGHT / 2)
    # mw = Window(x, y, DEFAULT_WIDTH, DEFAULT_HEIGHT, operating_system, cfg=cfg)
    mw = ControllingWindow(x, y, DEFAULT_WIDTH, DEFAULT_HEIGHT, operating_system, cfg=cfg, debug=True)
    mw.show()

    sys.exit(app.exec())
