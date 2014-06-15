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
    self.panel = wx.Panel(self, wx.ID_ANY)
    #img_path = os.path.abspath("../../../../../../53977.png")
    img_path = os.path.abspath("../../../../../../Pictures/duccati.png")
    bitmap = wx.Bitmap(img_path, type=wx.BITMAP_TYPE_PNG)
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)

    topSizer        = wx.BoxSizer(wx.VERTICAL)
    btnSizer        = wx.BoxSizer(wx.HORIZONTAL)
    
    okBtn = wx.Button(self.panel, wx.ID_ANY, 'OK')
    cancelBtn = wx.Button(self.panel, wx.ID_ANY, 'Cancel')
    self.Bind(wx.EVT_BUTTON, self.onOK, okBtn)
    self.Bind(wx.EVT_BUTTON, self.onCancel, cancelBtn)

    btnSizer.Add(okBtn, 0, wx.ALL, 5)
    btnSizer.Add(cancelBtn, 0, wx.ALL, 5)
    topSizer.Add(btnSizer, 0, wx.ALL|wx.CENTER, 5)
    self.panel.SetSizer(topSizer)
    topSizer.Fit(self)

  def onOK(self, event):
    # Do something
    print 'onOK handler'

  def onCancel(self, event):
    print "cancel clicked"
    #self.closeProgram()



if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()
