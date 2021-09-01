import json
from pathlib import Path
import os


class Configurator:
    DEBUG = None
    homepath = None
    cfg_file = None
    keys = {"contigsequences": "", "contigcoverages": "", "prodigal": "", "fetchMG": "", "fetchMGres": "",
            "kmeredata": ""}
    PRODIGAL_KEY = "prodigal"
    FETCHMG_KEY = "fetchMG"
    FETCHMGRES_KEY = "fetchMGres"
    CONTIGSEQ_KEY = "contigsequences"
    CONTIGCOV_KEY = "contigcoverages"
    KMERE_KEY = "kmeredata"

    def __init__(self, path, debug=False):
        self.homepath = path
        self.DEBUG = debug
        filepath = f"{path}{os.path.sep}paths.json"
        if os.path.exists(filepath):
            self.cfg_file = open(filepath, "r+")
            self.keys = json.load(self.cfg_file)
            self.is_new = False
            if self.DEBUG:
                print(f"[DEBUG] Configurator.__init__(): File exists.")
        else:
            self.cfg_file = open(filepath, "w+")
            self.cfg_file.write(json.dumps(self.keys))
            self.cfg_file.flush()
            self.is_new = True
            if self.DEBUG:
                print(f"[DEBUG] Configurator.__init__(): Created new File.")

    def write(self, key: str, val):
        self.keys[key] = val
        self.cfg_file.seek(0)
        self.cfg_file.truncate()
        self.cfg_file.write(json.dumps(self.keys))
        self.cfg_file.flush()

    def read(self, key: str):
        return self.keys[key]

    def on_close(self):
        if self.DEBUG:
            print("[DEBUG] Configurator.on_close()")
        self.cfg_file.close()
