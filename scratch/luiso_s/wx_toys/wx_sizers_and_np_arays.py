"""
Simple app that demonstrates how to use a wx.StaticBitmap, specifically
replacing bitmap dynamically.

Note: there needs to be an "Images" directory with one or more jpegs in it in the
      current working directory for this to work

Test most recently on OS-X wxPython 2.9.3.1

But it been reported to work on lots of other platforms/versions

"""
import wx, os, numpy
import matplotlib.pyplot as plt
class TestFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)

        # there needs to be an "Images" directory with one or more jpegs 
        # in it in the current working directory for this to work
        self.jpgs = GetJpgList("/home/lui/Pictures/") # get all the jpegs
                                                      # in the Images directory
        self.CurrentJpg = 0

        self.MaxImageSize = 500

        b = wx.Button(self, -1, "Next Reflection ")
        b_a = wx.Button(self, -1, "Previous Reflection")
        b.Bind(wx.EVT_BUTTON, self.DisplayNext)
        b_a.Bind(wx.EVT_BUTTON, self.DisplayPrev)

        # starting with an EmptyBitmap, the real one will get put there
        # by the call to .DisplayNext()
        self.Image = wx.StaticBitmap(self, bitmap=wx.EmptyBitmap(
                                     self.MaxImageSize, self.MaxImageSize
                                     ))

        self.DisplayNext()

        # Using a Sizer to handle the layout: I never  use absolute positioning
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(b, 0, wx.CENTER | wx.ALL,10)


        # adding stretchable space before and after centers the image.
        box.Add((1,1),1)
        box.Add(self.Image
                , 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL | wx.ADJUST_MINSIZE
                , 10)

        box.Add((1,1),1)
        box.Add(b_a, 0, wx.CENTER | wx.ALL,10)

        self.SetSizerAndFit(box)

        wx.EVT_CLOSE(self, self.OnCloseWindow)

    def DisplayNext(self, event = None):
        # load the image
        Img = wx.Image(self.jpgs[self.CurrentJpg], wx.BITMAP_TYPE_JPEG)

        # scale the image, preserving the aspect ratio
        W = Img.GetWidth()
        H = Img.GetHeight()
        if W > H:
            NewW = self.MaxImageSize
            NewH = self.MaxImageSize * H / W
        else:
            NewH = self.MaxImageSize
            NewW = self.MaxImageSize * W / H
        Img = Img.Scale(NewW,NewH)

        # convert it to a wx.Bitmap, and put it on the wx.StaticBitmap
        self.Image.SetBitmap(wx.BitmapFromImage(Img))

        # You can fit the frame to the image, if you want.
        self.Fit()
        self.Layout()
        self.Refresh()

        self.CurrentJpg += 1
        if self.CurrentJpg > len(self.jpgs) -1:
            self.CurrentJpg = 0

    def DisplayPrev(self, event = None):
        print "test"
        np_img = build_np_img(width = 164, height = 64)

        Img = GetBitmap_from_np_array(np_img)

        self.Image.SetBitmap(Img)
        self.Fit()
        self.Layout()
        self.Refresh()

    def OnCloseWindow(self, event):
        self.Destroy()

def GetBitmap_from_np_array(np_img_2d):
  fig = plt.figure()
  # remember to make sure this is our convention in (x, y) vs (row, col)
  plt.imshow(numpy.transpose(np_img_2d), interpolation = "nearest")
  fig.canvas.draw()
  width, height = fig.canvas.get_width_height()
  np_buf = numpy.fromstring ( fig.canvas.tostring_rgb(), dtype=numpy.uint8 )
  np_buf.shape = (width, height, 3)
  np_buf = numpy.roll(np_buf, 3, axis = 2)
  image = wx.EmptyImage(width, height)
  image.SetData( np_buf )
  #image.SetData( np_buf.tostring()) # looks like there is no need to convert
  wxBitmap = image.ConvertToBitmap()

  return wxBitmap

  
  
def build_np_img(width = 64, height = 64):
  data2d = numpy.zeros( (width, height), 'float')
  print "width, height =", width, height
  for x in range(0, width):
    for y in range(0, height):
      #data2d[x,y] = numpy.sqrt(x*x + y*y)
      data2d[x,y] = x + y
  data2d[width/4:width*3/4,height/4:height*3/4] = 0
  print "data2d.max =", data2d.max()
  return data2d
  
  
  
def GetJpgList(dir):
    jpgs = [f for f in os.listdir(dir) if f[-4:] == ".jpg"]
    # print "JPGS are:", jpgs
    return [os.path.join(dir, f) for f in jpgs]

class App(wx.App):
    def OnInit(self):

        frame = TestFrame(None, -1, "wxBitmap Test", wx.DefaultPosition
        ,(550,200))
        self.SetTopWindow(frame)
        frame.Show(True)
        return True

if __name__ == "__main__":
    app = App(0)
    app.MainLoop()