#!/usr/bin/python

import sys

import PySide
from PySide.QtGui import QApplication
from PySide.QtGui import QMessageBox

# Create the application object
app = QApplication(sys.argv)

# Create a simple dialog box
msgBox = QMessageBox()
msgBox.setText("Hello World - using PySide")
msgBox.exec_()