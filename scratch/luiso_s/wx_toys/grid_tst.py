# exercise with grids copied from:
# http://www.blog.pythonlibrary.org/2010/03/18/wxpython-an-introduction-to-grids/
# then evolved to our needs and to our coding style

import wx
import wx.grid as gridlib

########################################################################
class MyForm(wx.Frame):

  def __init__(self):
    """Constructor"""
    super(MyForm, self).__init__(parent=None, title="A Simple toying Grid")
    panel = wx.Panel(self)

    myGrid = gridlib.Grid(panel)
    myGrid.CreateGrid(12, 8)
    myGrid.SetCellValue(5, 0, "123")
    myGrid.EnableEditing(False)

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(myGrid, 1, wx.EXPAND)
    panel.SetSizer(sizer)

if __name__ == "__main__":
  app = wx.PySimpleApp()
  frame = MyForm().Show()
  app.MainLoop()

