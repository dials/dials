import numpy
import matplotlib

matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import wx

class RandomPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent)
    self._figure = Figure()
    self._axes = self._figure.add_subplot(1,1,1)
    self._canvas = FigureCanvas(self, -1, self._figure)
    self._sizer = wx.BoxSizer(wx.VERTICAL)
    self._sizer.Add(self._canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
    self.SetSizer(self._sizer)
    self.Fit()
    self._scale = 2.0
    self._canvas.Bind(wx.EVT_LEFT_DOWN, self.up)
    self._canvas.Bind(wx.EVT_RIGHT_DOWN, self.down)
    return

  def draw(self):
    size = 200

    im = numpy.empty([size, size])

    import random
    for i in range(size):
      for j in range(size):
        im[i,j] = random.random() * self._scale

    plot = self._axes.imshow(im)
    plot.set_clim(0.0, 1.0)
    self._canvas.draw()

    return

  def set_scale(self, scale):
    self._scale = scale
    self.draw()
    return

  def up(self, event):
    self._scale += 0.1
    self.draw()
    return

  def down(self, event):
    self._scale -= 0.1
    self.draw()
    return

if __name__ == "__main__":
  app = wx.PySimpleApp()
  frame = wx.Frame(None, title='Random Numbers')
  panel = RandomPanel(frame)
  panel.draw()
  frame.Show()
  app.MainLoop()
