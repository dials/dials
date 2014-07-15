import os
import wx
class MyApp(wx.App):
  def OnInit(self):
    self.frame = MyFrame(None, title="Bitmaps")
    self.SetTopWindow(self.frame)
    self.frame.Show()
    return True
class MyFrame(wx.Frame):
  def __init__(self, parent, id=wx.ID_ANY, title="",
      pos=wx.DefaultPosition, size=wx.DefaultSize,
      style=wx.DEFAULT_FRAME_STYLE,
      name="MyFrame"):
    super(MyFrame, self).__init__(parent, id, title,
                                  pos, size, style, name)

    self.text01 = wx.TextCtrl(self, -1, "enter number", size = (300, -1))
    btn_read = wx.Button(self, -1, "Read number")

    h_box = wx.BoxSizer(wx.HORIZONTAL)
    v_box = wx.BoxSizer(wx.VERTICAL)
    v_box.Add(self.text01, 0, wx.CENTER | wx.ALL,5)
    v_box.Add(btn_read, 0, wx.CENTER | wx.ALL,5)
    h_box.Add(v_box)
    btn_read.Bind(wx.EVT_BUTTON, self.read_num_from_txt)

    self.SetSizerAndFit(h_box)

  def read_num_from_txt(self, event = None):
    print "clicked"
    a = int(self.text01.GetValue())
    print "a =", a
    print "a * 2 =", a * 2
if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()
