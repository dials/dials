import wx

class RdButtonFrm(wx.Frame):
  def __init__(self):
    wx.Frame.__init__(self, None, -1, 'Radio Example', size = (200, 200))
    panel = wx.Panel(self, -1)
    radio1= wx.RadioButton(panel, -1, "Elmo", pos= (20, 50), style = wx.RB_GROUP)
    radio2= wx.RadioButton(panel, -1, "Ern", pos= (20, 80))#, style = wx.RB_GROUP)
    radio3= wx.RadioButton(panel, -1, "Brt", pos= (20, 110))#, style = wx.RB_GROUP)
    text1 = wx.TextCtrl(panel, -1, "", pos = (80, 50))
    text2 = wx.TextCtrl(panel, -1, "", pos = (80, 80))
    text3 = wx.TextCtrl(panel, -1, "", pos = (80, 110))
    self.texts = {"Elmo": text1, "Ern": text2, "Brt": text3}
    for eachText in [text2, text3]:
      eachText.Enable(False)
    for eachRadio in [radio1, radio2, radio3]:
      self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio, eachRadio)
    self.selectedText = text1

  def OnRadio(self, event):
    if self.selectedText:
      self.selectedText.Enable(False)
    radioSelected = event.GetEventObject()
    text = self.texts[radioSelected.GetLabel()]
    text.Enable(True)
    self.selectedText = text
  def tst(self, a):
    print "from tst"
    print a
    return 2*a

if __name__ == '__main__':
  app = wx.PySimpleApp()
  frst = RdButtonFrm()
  frst.Show()
  scnd = RdButtonFrm()
  scnd.Show()
  dos_veses = frst.tst(5)
  app.MainLoop()
