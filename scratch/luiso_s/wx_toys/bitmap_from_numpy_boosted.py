import wx
import numpy as np

def GetBitmap(width = 255, height = 32 ):
        array = np.zeros( (height, width, 3),'uint8')


        array[:,:,:] = [ (itst, itst, itst) for itst in range(width) ]


        image = wx.EmptyImage(width,height)
        image.SetData( array.tostring())
        wxBitmap = image.ConvertToBitmap()       # OR:  wx.BitmapFromImage(image)
        return wxBitmap


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

    self.panel = wx.Panel(self)
    bitmap = GetBitmap()
    self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if(__name__ == "__main__"):
  print ("Ini App ")
  app = MyApp(redirect=False)
  print ("running App")
  app.MainLoop()
