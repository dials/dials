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
    # Attributes

    text01 = wx.TextCtrl(self, -1, "enter txt", size = (300, -1))
    btn_nxt_refl = wx.Button(self, -1, "Previous Reflection")

    h_box = wx.BoxSizer(wx.HORIZONTAL)
    v_box = wx.BoxSizer(wx.VERTICAL)
    v_box.Add(text01, 0, wx.CENTER | wx.ALL,5)
    v_box.Add(btn_nxt_refl, 0, wx.CENTER | wx.ALL,5)
    h_box.Add(v_box)

    self.SetSizerAndFit(h_box)
    btn_nxt_refl.Bind(wx.EVT_BUTTON, self.DisplayNext_refl)
  def DisplayNext_refl(self, event = None):
    print "clicked"
    #'GetValue'
if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()
