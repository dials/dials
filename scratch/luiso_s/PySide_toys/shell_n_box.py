
from subprocess import call as shell_func
import sys
from PySide import QtGui


class Example(QtGui.QWidget):

    def __init__(self):
        super(Example, self).__init__()
        self.initUI()

    def initUI(self):

        self.btn = QtGui.QPushButton('Go Do it', self)
        self.btn.move(20, 20)
        self.btn.clicked.connect(self.clickeado)

        self.lin_txt = QtGui.QLineEdit(self)
        self.lin_txt.move(130, 22)

        self.setGeometry(300, 300, 290, 150)
        self.setWindowTitle('Shell dialog')
        self.show()

    def clickeado(self):

        shell_str = str(self.lin_txt.text())
        shell_func(shell_str, shell=True)

        self.lin_txt.setText(str(""))


def main():

    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
