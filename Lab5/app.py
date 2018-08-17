import sys
import csv
import os
import itertools as it
import functools as ft
from typing import *
from os.path import join
from PyQt5.QtWidgets import QDialog, QApplication
from guillya import Ui_MainWindow

class Aggregator(NamedTuple):
    name:        str
    connections: int
    price:       float

class Cloud(NamedTuple):
    name:  str
    price: int

class Connection(NamedTuple):
    name:         str
    refresh_time: float
    price:        float

class MeteoModel(NamedTuple):
    name:      str
    precision: int
    farthest:  int
    time:      int

class Meteostation(NamedTuple):
    name:  str
    temp:  float
    wind:  float
    hum: float
    press: float
    price: float

class Storage(NamedTuple):
    name:     str
    max_time: int
    price:    int

class Configuration(NamedTuple):
    aggregator: Aggregator
    cloud:      Cloud
    connection: Connection
    model:      MeteoModel
    station:    Meteostation
    storage:    Storage

    def __str__(self):
        return  ( "\t".join(map(lambda x: x.name, self))
                + "\t" + str(self.price)
                + "\t" + str(self.precision)
                )
    @property
    def precision(self):
        s = self
        return ( s.aggregator.connections
               * s.connection.refresh_time
               * ( s.station.temp
                 + s.station.wind
                 + s.station.hum
                 + s.station.press
                 )
               * s.model.precision
               * s.model.farthest
               * s.storage.max_time
               )

    @property
    def price(self):
        s = self
        return ( s.aggregator.price
               + s.cloud.price
               + s.connection.price
               + s.station.price
               + s.storage.price
               )

def get_variants(path: str, fls: Dict[NamedTuple, str]):
    res = {k: [] for k in fls.keys()}
    for k, v in fls.items():
        locres = []
        with open(resource_path(join(path, v)), 'r', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            next(reader) # skipping header
            for row in reader:
                if   k == Aggregator:
                    n, c, p = row
                    locres.append(k(name        = n,
                                    connections = int(c),
                                    price       = float(p),
                                   )
                                 )
                elif k == Cloud:
                    n, p = row
                    locres.append(k(name  = n,
                                    price = float(p),
                                   )
                                 )
                elif k == Connection:
                    n, rt, p = row
                    locres.append(k(name         = n,
                                    refresh_time = float(rt),
                                    price        = float(p),
                                   )
                                 )
                elif k == MeteoModel:
                    n, pr, f, t = row
                    locres.append(k(name      = n,
                                    precision = int(pr),
                                    farthest  = int(f),
                                    time      = int(t),
                                   )
                                 )
                elif k == Meteostation:
                    n, t, w, h, pr, p = row
                    locres.append(k(name  = n,
                                    temp  = float(t),
                                    wind  = float(w),
                                    hum   = float(h),
                                    press = float(pr),
                                    price = float(p)
                                   )
                                 )
                elif k == Storage:
                    n, t, p = row
                    locres.append(k(name      = n,
                                    max_time  = int(t),
                                    price     = int(p),
                                   )
                                 )
        res[k] = locres
    return res

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return join(base_path, relative_path)

class AppWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.show()  
        self.ui.pushButton.clicked.connect(self.DoTheStuff)
        parts_path = "parts"
        files = {
            Aggregator:   "aggr.csv",
            Cloud:        "cloud-platform.csv",
            Connection:   "conn.csv",
            MeteoModel:   "meteomodel.csv",
            Meteostation: "meteostation.csv",
            Storage:      "storage.csv",
        }
        self.variants = get_variants(parts_path, files)

    def DoTheStuff(self):
        agreg_con = (self.ui.lineEdit.text())

        speed_renew = (self.ui.lineEdit_2.text())

        acc_temp = (self.ui.lineEdit_3.text())
        acc_wind = (self.ui.lineEdit_4.text())
        acc_hum = (self.ui.lineEdit_5.text())
        acc_pres = (self.ui.lineEdit_6.text())

        acc_forec = (self.ui.lineEdit_7.text())

        time_forec = (self.ui.lineEdit_8.text())
        time_to_store = (self.ui.lineEdit_9.text())
        provider  = (self.ui.comboBox.currentIndex())
        provider = self.variants[Cloud][provider]

        # low_price = float(self.ui.lineEdit_10.text())
        high_price = float(self.ui.lineEdit_11.text())

        low_accur = float(self.ui.lineEdit_12.text())
        # high_accur = float(self.ui.lineEdit_13.text())
        all_combinations = it.starmap(Configuration, it.product(*self.variants.values()))
        ac = lambda c: (   c.aggregator.connections >= int(agreg_con)
                       and c.connection.refresh_time <= float(speed_renew)
                       and c.station.temp <= float(acc_temp)
                       and c.station.wind <= float(acc_wind)
                       and c.station.hum <= float(acc_hum)
                       and c.station.press <= float(acc_pres)
                       and c.model.precision <= float(acc_forec)
                       and c.model.farthest >= int(time_forec)
                       and c.storage.max_time >= int(time_to_store)
                       and c.cloud == provider
                       )
        fc = lambda c: (   low_accur <= c.precision
                       and c.price     <= high_price
                       )
        oldres = list(filter(lambda c: fc(c) and ac(c), all_combinations))
        # pareto = "\n※".join(map(str, res))
        seen = Configuration(*(set() for i in oldres[0])) if oldres else None
        in_sets = lambda config, sets: any((c in s) for (c,s) in zip(config, sets))
        add_in_sets = lambda config, sets: [s.add(c) for (c,s) in zip(config, sets)]
        res = [c for c in oldres if not in_sets(c, seen) and add_in_sets(c, seen)]
        parset = Configuration(*(set() for i in oldres[0])) if oldres else None
        pareto = list(map(lambda config: add_in_sets(config, parset),
                          oldres))
        self.ui.plainTextEdit.setPlainText( (("#Pareto set:\n"
                                            + "※"
                                            + "\n※ ".join(map(lambda x: str(",".join(map(str, x))),
                                                              parset))

                                            + "\n"
                                            ) if seen else "")
                                           + (("#Results:\n"
                                             + "※" + "\n※".join(map(str, res))
                                             + "\n") if res else "No results, sorry!")
                                          )

app = QApplication(sys.argv)
w = AppWindow()
w.show()
sys.exit(app.exec_())
