import json
from pathlib import Path
import os


class Configurator:
    DEBUG = True
    homepath = None
    cfg_file = None
    keys = {"prodigal": "", "fetchMG": "", "working": ""}
    PRODIGAL_KEY = "prodigal"
    FETCHMG_KEY = "fetchMG"
    DATA_KEY = "working"

    def __init__(self, path):
        self.homepath = path
        filepath = f"{path}/paths.json"
        if os.path.exists(filepath):
            self.cfg_file = open(filepath, "r+")
            self.keys = json.load(self.cfg_file)
            print("existing")
        else:
            self.cfg_file = open(filepath, "w+")
            self.cfg_file.write(json.dumps(self.keys))
            self.cfg_file.flush()
            print("written")

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
