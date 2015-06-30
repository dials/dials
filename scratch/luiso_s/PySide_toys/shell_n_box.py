from subprocess import call as shell_func
import sys
from PySide import QtGui

class Example(QtGui.QWidget):

    def __init__(self):
        super(Example, self).__init__()
        self.initUI()

    def initUI(self):

        self.btn1 = QtGui.QPushButton('    import    ', self)
        self.btn1.move(20, 20)
        self.btn1.clicked.connect(self.B_clicked1)
        self.btn2 = QtGui.QPushButton(' find_spots ', self)
        self.btn2.move(20, 60)
        self.btn2.clicked.connect(self.B_clicked2)
        self.btn3 = QtGui.QPushButton('    index     ', self)
        self.btn3.move(20, 100)
        self.btn3.clicked.connect(self.B_clicked3)

        self.btn_go = QtGui.QPushButton('      Go      ', self)
        self.btn_go.move(300,180)
        self.btn_go.clicked.connect(self.B_go_clicked)

        self.lin_txt = QtGui.QLineEdit(self)
        self.lin_txt.move(20, 182)

        self.setGeometry(100, 200, 500, 250)
        self.setWindowTitle('Shell dialog')
        self.show()

    def B_clicked1(self):
        self.lin_txt.setText(str("dials.import"))

    def B_clicked2(self):
        self.lin_txt.setText(str("dials.find_spots"))

    def B_clicked3(self):
        self.lin_txt.setText(str("dials.index"))

    def B_go_clicked(self):
        shell_str = str(self.lin_txt.text())
        shell_func(shell_str, shell=True)
        self.lin_txt.setText(str(""))

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
