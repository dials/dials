from __future__ import division

class spot_wrapper(object):
  def __init__(self, params):
    self.params = params
    self.frame = None

  def display(self, imagesets, reflections, crystals=None):
    import wx
    from dials.util.spotfinder_frame import SpotFrame

    app = wx.App()

    self.frame = SpotFrame(None, -1, "X-ray image display", size=(800,720),
      pos=(100,100),
      params=self.params,
      imagesets=imagesets,
      reflections=reflections,
      crystals=crystals)
    self.frame.SetSize((1024,780))
    self.frame.Show()

    from rstbx.slip_viewer.frame import chooser_wrapper

    for imageset in imagesets:
      for idx in xrange(len(imageset.indices())):
        self.frame.add_file_name_or_data(chooser_wrapper(imageset, idx))
    self.frame.load_image(chooser_wrapper(imagesets[0], 0))

    app.MainLoop()

  def load_image(self, filename):
    from dials.util.spotfinder_frame import create_load_image_event
    if self.frame is not None:
      create_load_image_event(self.frame, filename)
