"""
@author: jez2bp
"""


import _converter.conversion.pltparser
from pathlib import Path
import json
import os

class ReadPLT:
    """The class provides functions and methods for passing between pltparser.py
    and convert_main.py
    """

    def __init__(self, pltpath):
        self.pltpath = pltpath

    def extract_plt_data(self):
        """The function give back the necessary data from the plt parseing
        """
        pltdict = {}
        sourcenames = []
        with open(os.path.join(Path(os.path.dirname(os.path.abspath(__file__))).parent.parent.parent, '_temp',
                               'config.json')) as f:
            config = json.load(f)
        for signals in config["signal_list"]:
            if len(signals) == 1:
                pltdict[signals[0]] = {"name": signals[0]}
                pltdict[signals[0]]["exp"] = None
                sourcenames.append(signals[0])
            elif len(signals) == 2:
                pltdict[signals[-1]] = {"name": signals[0]}
                pltdict[signals[-1]]["exp"] = None
                sourcenames.append(signals[0])
        return sourcenames, pltdict
