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

    img_path = os.path.abspath("../../../../../../Pictures/dials_logo01.png")
    #img_path = os.path.abspath("../../../../../../Pictures/duccati.png")

    self.panel = wx.Panel(self, wx.ID_ANY)

    bitmap = wx.Bitmap(img_path, type=wx.BITMAP_TYPE_PNG)
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)

    self.build_content()
    self.resize()

  def build_content(self):
    topSizer        = wx.BoxSizer(wx.VERTICAL)
    btnSizer        = wx.BoxSizer(wx.HORIZONTAL)

    prev_Btn = wx.Button(self.panel, wx.ID_ANY, 'Prev Refl')
    next_Btn = wx.Button(self.panel, wx.ID_ANY, 'Next Refl')
    self.Bind(wx.EVT_BUTTON, self.on_prev, prev_Btn)
    self.Bind(wx.EVT_BUTTON, self.on_next, next_Btn)
    topSizer.Add(self.bitmap, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 3)
    #box.Add(self.Image, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE, 10)
    btnSizer.Add(prev_Btn, 0, wx.ALL, 3)
    btnSizer.Add(next_Btn, 0, wx.ALL, 3)
    topSizer.Add(btnSizer, 0, wx.ALL|wx.CENTER, 3)
    self.my_main_sizer = topSizer

  def resize(self):

    self.panel.SetSizer(self.my_main_sizer)


    #self.panel.SetAutoLayout(True)
    self.panel.Refresh()
    self.panel.Update()
    self.panel.Layout()
    self.Update()
    self.Layout()

    self.my_main_sizer.Fit(self)




  def on_prev(self, event):
    print 'on_prev handler'
    img_path = os.path.abspath("../../../../../../Pictures/dials_logo02.png")
    new_bitmap = wx.Bitmap(img_path, type=wx.BITMAP_TYPE_PNG)
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=new_bitmap)
    self.resize()

  def on_next(self, event):
    print "on_next handler"
    img_path = os.path.abspath("../../../../../../Pictures/dials_logo03.png")
    new_bitmap = wx.Bitmap(img_path, type=wx.BITMAP_TYPE_PNG)
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=new_bitmap)
    self.resize()



if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()
