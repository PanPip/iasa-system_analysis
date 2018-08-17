import sys
from PyQt5.QtWidgets import QDialog, QApplication, QTableWidgetItem
from guillya6 import Ui_MainWindow
from mpart import calc

class AppWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.show()  
        self.ui.pushButton.clicked.connect(self.PrintResults)
        self.ui.pushButton_2.clicked.connect(lambda : self.ui.stackedWidget.setCurrentIndex(0))
        self.ui.pushButton_3.clicked.connect(lambda : self.ui.stackedWidget.setCurrentIndex(1))
        self.ui.pushButton_4.clicked.connect(lambda : self.ui.stackedWidget.setCurrentIndex(2))

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
        low_price = (self.ui.lineEdit_10.text())
        high_price = (self.ui.lineEdit_11.text())
        low_accur = (self.ui.lineEdit_12.text())
        high_accur = (self.ui.lineEdit_13.text())
        self.ui.plainTextEdit.setPlainText(str(provider))

    def PrintResults(self):
        res = calc()
        self.ui.tableWidget_6.setItem(0,1, QTableWidgetItem(str(round((res[0]),3))))
        self.ui.tableWidget_6.setItem(1,1, QTableWidgetItem(str(round((res[1]),3)))) 
        self.ui.tableWidget_6.setItem(2,1, QTableWidgetItem(str(round((res[2]),3))))
        self.ui.tableWidget_6.setItem(3,1, QTableWidgetItem(str(round((res[3]),3))))
        self.ui.tableWidget_6.setItem(4,1, QTableWidgetItem(str(round((res[4]),3))))
        self.ui.tableWidget_6.setItem(5,1, QTableWidgetItem(str(round((res[5]),3)))) 

app = QApplication(sys.argv)
w = AppWindow()
w.show()
sys.exit(app.exec_())
