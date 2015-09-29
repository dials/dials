from __future__ import division
class spot_wrapper(object):
  def __init__(self, params):
    self.params = params

  def display(self, imagesets, reflections, crystals=None):
    import wx
    from dials.util.spotfinder_frame import SpotFrame

    app = wx.App(0)

    frame = SpotFrame(None, -1, "X-ray image display", size=(800,720),
      pos=(100,100),
      params=self.params,
      imagesets=imagesets,
      reflections=reflections,
      crystals=crystals)
    frame.SetSize((1024,780))
    frame.Show()

    from rstbx.slip_viewer.frame import chooser_wrapper

    i = 0
    for imageset in imagesets:
      for idx in imageset.indices():
        i += 1
        #print i
        frame.add_file_name_or_data(chooser_wrapper(imageset, idx))
    idx = imagesets[0].indices()[0]
    frame.load_image(chooser_wrapper(imagesets[0],idx))

    app.MainLoop()
