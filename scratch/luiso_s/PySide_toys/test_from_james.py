
from __future__ import division
from PyQt4.QtGui import QMainWindow
from PyQt4.QtGui import QApplication
from PyQt4.QtGui import QTreeView
from PyQt4.QtGui import QVBoxLayout
from PyQt4.QtGui import QPushButton
from PyQt4.QtGui import QWidget
from PyQt4.QtGui import QStandardItemModel
from PyQt4.QtGui import QSortFilterProxyModel
from PyQt4.QtGui import QStandardItem
from PyQt4.QtGui import QAbstractItemView
from PyQt4.QtGui import QHeaderView
from PyQt4.QtGui import QStyledItemDelegate
from PyQt4.QtGui import QSpinBox
from PyQt4.QtGui import QDoubleSpinBox
from PyQt4.QtGui import QLineEdit
from PyQt4.QtGui import QComboBox
from PyQt4.QtCore import Qt
from PyQt4.QtCore import QSize
from PyQt4.QtCore import QString
from PyQt4.QtCore import QRegExp


class StringEditor(QLineEdit):

  def __init__(self, parent=None):
    super(StringEditor, self).__init__(parent)

  def setEditorData(self, index):
    value = index.model().data(index, Qt.EditRole).toString()
    self.setText(value)
    
  def setModelData(self, model, index):
    value = self.text()
    model.setData(index, value, Qt.EditRole)


class IntEditor(QSpinBox):

  def __init__(self, parent=None, value_min=None, value_max=None):
    super(IntEditor, self).__init__(parent)
    if value_min is not None:
      self.setMinimum(value_min)
    else:
      self.setMinimum(-1e9)
    if value_max is not None:
      self.setMaximum(value_max)
    else:
      self.setMaximum(1e9)

  def setEditorData(self, index):
    value = index.model().data(index, Qt.EditRole).toInt()
    self.setValue(value[0])

  def setModelData(self, model, index):
    self.interpretText()
    value = self.value()
    model.setData(index, value, Qt.EditRole)


class FloatEditor(QDoubleSpinBox):

  def __init__(self, parent=None, value_min=None, value_max=None):
    super(FloatEditor, self).__init__(parent)
    if value_min is not None:
      self.setMinimum(value_min)
    else:
      self.setMinimum(-1e9)
    if value_max is not None:
      self.setMaximum(value_max)
    else:
      self.setMaximum(1e9)
    self.setDecimals(2)
  
  def setEditorData(self, index):
    value = index.model().data(index, Qt.EditRole).toFloat()
    self.setValue(value[0])
   
  def setModelData(self, model, index):
    self.interpretText()
    value = self.value()
    model.setData(index, value, Qt.EditRole)


class ChoiceEditor(QComboBox):

  def __init__(self, parent=None, choices=None):
    super(ChoiceEditor, self).__init__(parent)
    if choices is not None:
      for choice in choices:
        self.addItem(choice)
  
  def setEditorData(self, index):
    value = index.model().data(index, Qt.EditRole).toString()
    self.setCurrentIndex(self.findText(value))
   
  def setModelData(self, model, index):
    value = self.currentText()
    model.setData(index, value, Qt.EditRole)


class BoolEditor(QComboBox):

  def __init__(self, parent=None):
    super(BoolEditor, self).__init__(parent)
    self.addItem("False")
    self.addItem("True")

  def setEditorData(self, index):
    value = index.model().data(index, Qt.EditRole).toString()
    self.setCurrentIndex(self.findText(value))

  def setModelData(self, model, index):
    value = self.currentText()
    model.setData(index, value, Qt.EditRole)


class ParameterWindowDelegate(QStyledItemDelegate):

  def __init__(self, parent=None):

    super(ParameterWindowDelegate, self).__init__(parent)

  def createEditor(self, parent, option, index):
    parameter = index.model().data(index, Qt.UserRole+1).toPyObject()
    dtype = parameter.type
    ptype = dtype.phil_type
    if ptype == 'str':
      editor = StringEditor(parent)
    elif ptype == 'float':
      editor = FloatEditor(parent, dtype.value_min, dtype.value_max)
    elif ptype == 'int':
      editor = IntEditor(parent, dtype.value_min, dtype.value_max)
    elif ptype == 'choice':
      def strip(w):
        w = str(w)
        if w.startswith("*"):
          return w[1:]
        return w
      choices = [strip(w) for w in parameter.words]
      editor = ChoiceEditor(parent, choices)
    elif ptype == 'bool':
      editor = BoolEditor(parent)
    else:
      raise RuntimeError("Handle type %s" % dtype)
    return editor

  def setEditorData(self, editor, index):
    editor.setEditorData(index)

  def setModelData(self, editor, model, index):
    editor.setModelData(model, index)

  def updateEditorGeometry(self, editor, option, index):
    editor.setGeometry(option.rect)

  def sizeHint(self, option, index):
    size = super(ParameterWindowDelegate, self).sizeHint(option, index)
    size.setHeight(size.height() * 2)
    return size


class ParameterTreeModel(QStandardItemModel):

  def __init__(self):
    super(ParameterTreeModel, self).__init__()

  def data(self, index, role=Qt.DisplayRole):
    if role == Qt.ToolTipRole:
      parameter = index.model().data(index, Qt.UserRole+1).toPyObject()
      if parameter is None:
        return ""
      return str(parameter.help)
    return super(ParameterTreeModel, self).data(index, role)


class ParameterWindow(QWidget):

  def __init__(self, parent=None, parameters=None, expert_level=0):
    
    # Init the parent 
    super(ParameterWindow, self).__init__(parent)

    # Save the parameters
    self.parameters = parameters

    # Create the model
    model = ParameterTreeModel()
    model.setHorizontalHeaderLabels(['Parameter', 'Value'])

    # Traverse the parameter tree and add to the tree view widget
    # FIXME Store data in tree nodes better
    def visit(root, parameter):
      if parameter.is_scope:
        name_node = QStandardItem(parameter.name)
        name_node.setFlags(Qt.NoItemFlags | Qt.ItemIsEnabled)
        name_node.setData(parameter, Qt.UserRole + 1)
        for obj in parameter.objects:
          visit(name_node, obj)
        root.appendRow(name_node)
      elif parameter.is_definition:
        name_node = QStandardItem(parameter.name)
        name_node.setFlags(Qt.NoItemFlags | Qt.ItemIsEnabled)
        name_node.setData(parameter, Qt.UserRole + 1)
        value_node = QStandardItem(str(parameter.extract()))
        value_node.setData(parameter, Qt.UserRole + 1)
        root.appendRow([name_node, value_node])
      else:
        raise RuntimeError('Handle This!')

    # Populate the tree
    if parameters is not None:
      for obj in parameters.objects:
        visit(model, obj)

    # Create a parameter tree
    self.tree = QTreeView()
    self.tree.setItemDelegate(ParameterWindowDelegate())
    self.tree.setAlternatingRowColors(True)
    self.tree.setSortingEnabled(False)
    self.tree.setHeaderHidden(False)
    self.tree.setSelectionBehavior(QAbstractItemView.SelectItems)

    # Set the model
    proxy_model = QSortFilterProxyModel()
    proxy_model.setSourceModel(model)
    proxy_model.setFilterKeyColumn(0)
    self.tree.setModel(proxy_model)
    self.tree.header().setStretchLastSection(True)
    self.tree.header().setResizeMode(0, QHeaderView.ResizeToContents)

    # Start everything expanded
    self.tree.expandAll()

    # Create the layout
    layout = QVBoxLayout()
    layout.setMargin(0)
    layout.addWidget(self.tree)
    self.setLayout(layout)

  def setExpertLevel(self, level):
    pass

  def setSearchString(self, text):
    print text
    regex = QRegExp(".*%s.*" % text, Qt.CaseInsensitive, QRegExp.FixedString)
    print regex
    self.tree.model().setFilterRegExp(regex)


class IntegrateParameterWindow(ParameterWindow):

  def __init__(self, parent=None, expert_level=0):
    from dials.command_line.integrate import phil_scope

    # Extract the parameters
    params = phil_scope.extract()

    # Init parent
    super(IntegrateParameterWindow, self).__init__(
      parent, 
      phil_scope,
      expert_level)

    
class MainWindow(QMainWindow):

  def __init__(self, parent=None):
   
    # Call the parent constructor
    super(MainWindow, self).__init__(parent)

    # Create the parameter window widget
    params = IntegrateParameterWindow(expert_level=0)

    # Create the search widget
    search = QLineEdit()
    search.setPlaceholderText("Search...")
    search.textChanged.connect(params.setSearchString)

    # Create the expert level widget
    expert = QSpinBox()
    expert.setMinimum(0)
    expert.setMaximum(1)
    expert.valueChanged.connect(params.setExpertLevel)

    # Create the window layout
    layout = QVBoxLayout()
    layout.addWidget(search)
    layout.addWidget(params)
    layout.addWidget(expert)

    # Setup the window contents 
    window = QWidget()
    window.setLayout(layout)
    self.setCentralWidget(window)


if __name__ == '__main__':
  import sys

  # Create the application
  app = QApplication(sys.argv)

  # Create the main window
  window = MainWindow()
  window.resize(800, 600)
  window.show()

  # Execute the application
  sys.exit(app.exec_())
