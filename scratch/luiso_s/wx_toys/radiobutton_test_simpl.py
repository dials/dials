import wx

class RdButtonFrm(wx.Frame):
  def __init__(self):
    wx.Frame.__init__(self, None, -1, 'Radio Example', size = (200, 200))
    panel = wx.Panel(self, -1)
    radio1= wx.RadioButton(panel, -1, "Elmo", pos= (20, 50), style = wx.RB_GROUP)
    radio2= wx.RadioButton(panel, -1, "Ern", pos= (20, 80))
    radio3= wx.RadioButton(panel, -1, "Brt", pos= (20, 110))

    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio1, radio1)
    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio2, radio2)
    self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio3, radio3)

  def OnRadio1(self, event = None):
    print "clicked on Radio1"
  def OnRadio2(self, event = None):
    print "clicked on Radio2"
  def OnRadio3(self, event = None):
    print "clicked on Radio3"

if __name__ == '__main__':
  app = wx.PySimpleApp()
  RdButtonFrm().Show()
  app.MainLoop()