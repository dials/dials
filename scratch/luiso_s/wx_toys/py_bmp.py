import os
import wx
class MyApp(wx.App):
    def OnInit(self, a_var = 1):
        self.frame = MyFrame(None, title="Bitmaps", b_var = a_var)
        self.SetTopWindow(self.frame)
        self.frame.Show()
        print "a_var =", a_var
        return True
class MyFrame(wx.Frame):
    def __init__(self, parent, id=wx.ID_ANY, title="",
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 style=wx.DEFAULT_FRAME_STYLE,
                 name="MyFrame", b_var = 0):
        super(MyFrame, self).__init__(parent, id, title,
                                      pos, size, style, name)
        # Attributes
        self.panel = wx.Panel(self)
        my_path = "./img_tmp.png"
        img_path = os.path.abspath(my_path)
        bitmap = wx.Bitmap(img_path, type=wx.BITMAP_TYPE_PNG)
        self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)
        # test var
        print "b_var =", b_var
if(__name__ == "__main__"):

    app = MyApp(False, 5)
    app.MainLoop()
