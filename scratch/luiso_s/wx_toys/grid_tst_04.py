# exercise with grids copied from:
# http://www.blog.pythonlibrary.org/2010/03/18/wxpython-an-introduction-to-grids/
# then evolved to our needs and to our coding style

import wx
import wx.grid as gridlib


class TupTable(gridlib.PyGridTableBase):
  def __init__(self, data, rowLabels=None, colLabels=None):
    gridlib.PyGridTableBase.__init__(self)
    self.data = data
    self.rowLabels = rowLabels
    self.colLabels = colLabels

  def GetNumberRows(self):
    return len(self.data)

  def GetNumberCols(self):
    return len(self.data[0])

  def GetColLabelValue(self, col):
    if self.colLabels:
      return self.colLabels[col]

  def GetRowLabelValue(self, row):
    if self.rowLabels:
      return self.rowLabels[row]

  def IsEmptyCell(self, row, col):
    return False

  def GetValue(self, row, col):
    return self.data[row][col]

  def SetValue(self, row, col, value):
    pass


class MyGrid(gridlib.Grid):
  def __init__(self, parent):
    """Constructor"""

    data = (("A", "B"),
            ("C", "D"),
            ("E", "Fxx123xx"),
            ("G", "G"),
            ("F", "F"),
            ("Q", "Q"))
    colLabels = ("Last", "Test")
    rowLabels = ("1", "2", "X", "40000", "E", "6")

    super(MyGrid, self).__init__(parent)
    tableBase = TupTable(data, rowLabels, colLabels)
    self.SetTable(tableBase)

    self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)


  def OnCellLeftClick(self, evt):
    print "OnCellLeftClick: (%d,%d) %s\n" % (evt.GetRow(),
                                             evt.GetCol(),
                                             evt.GetPosition())
    evt.Skip()

class MyForm(wx.Frame):

  def __init__(self):
    """Constructor"""
    super(MyForm, self).__init__(parent=None, title="An Grid for toying")
    panel = wx.Panel(self)
    #grid = SimpleGrid(self)
    myGrid = MyGrid(panel)

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(myGrid, 1, wx.EXPAND)
    panel.SetSizer(sizer)

if __name__ == "__main__":
  app = wx.PySimpleApp()
  frame = MyForm().Show()
  app.MainLoop()
