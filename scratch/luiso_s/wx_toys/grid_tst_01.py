# exercise with grids copied from:
# http://www.blog.pythonlibrary.org/2010/03/18/wxpython-an-introduction-to-grids/
# then evolved to our needs and to our coding style

import wx
import wx.grid as gridlib

class MyGrid(gridlib.Grid):

  def __init__(self, parent):
    """Constructor"""
    super(MyGrid, self).__init__(parent)
    self.CreateGrid(12, 8)
    self.SetColLabelValue(5, "test")
    # test all the events
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

    myGrid = MyGrid(panel)

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(myGrid, 1, wx.EXPAND)
    panel.SetSizer(sizer)

if __name__ == "__main__":
  app = wx.PySimpleApp()
  frame = MyForm().Show()
  app.MainLoop()
