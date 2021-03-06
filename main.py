import sys
import os
import platform
from lib import *
from qtgui import *
from PyQt5.QtWidgets import *

DEFAULT_WIDTH = 1280
DEFAULT_HEIGHT = 720
DEBUG = False
TEST = False

if __name__ == '__main__':
    app = QApplication(sys.argv)

    arguments = ['-d', '-t', '-h']  # -d: Debug True, -t: Test true
    ap = ArgParser(arguments)
    given_args = ap.parse(sys.argv)
    DEBUG = given_args['-d']
    TEST = given_args['-t']
    HELP = given_args['-h']

    if HELP:
        print("Possible Arguments: \n"
              "-d : DEBUG Mode\n"
              "-t : TEST Mode (Taking labels as additional input for testing performance of this software)\n"
              "-h : HELP (display this List)")
    else:
        operating_system = platform.system()

        working_path = None
        if operating_system == 'Linux':
            working_path = os.path.expanduser("~/BinMiner")
        elif operating_system == 'Windows':
            working_path = os.path.expandvars(R"C:\Users\$USERNAME\Documents\BinMiner")

        if not os.path.exists(working_path):
            os.makedirs(working_path)
        elif DEBUG:
            print(f"[DEBUG]{working_path} already exists.")
        if working_path:
            cfg = Configurator(working_path, DEBUG)
            size = app.primaryScreen().size()
            x = (size.width() / 2) - (DEFAULT_WIDTH / 2)
            y = (size.height() / 2) - (DEFAULT_HEIGHT / 2)
            mw = ControllingWindow(x, y, DEFAULT_WIDTH, DEFAULT_HEIGHT, cfg=cfg, debug=DEBUG, test=TEST)
            mw.show()

            sys.exit(app.exec())
        else:
            print(f"[ERROR] Unknown operating system or pathing problems.")
