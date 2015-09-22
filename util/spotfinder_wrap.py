from __future__ import division
class spot_wrapper(object):
  def __init__(self, params):
    self.params = params

  def display(self, imagesets, reflections, crystals=None):
    import wx
    from dials.util.spotfinder_frame import SpotFrame

    app = wx.App(0)

    # XXX Hacky workaround wxPython 3.0 crash on Linux
    # wx._core.PyAssertionError: C++ assertion "m_window" failed at ./src/gtk/dcclient.cpp(2043) in DoGetSize(): GetSize() doesn't work without window
    app.SetAssertMode(wx.PYAPP_ASSERT_SUPPRESS)

    frame = SpotFrame(None, -1, "X-ray image display", size=(800,720),
      pos=(100,100),
      params=self.params,
      imagesets=imagesets,
      reflections=reflections,
      crystals=crystals)
    frame.SetSize((1024,780))

    from rstbx.slip_viewer.frame import chooser_wrapper

    i = 0
    for imageset in imagesets:
      for idx in imageset.indices():
        i += 1
        #print i
        frame.add_file_name_or_data(chooser_wrapper(imageset, idx))
    idx = imagesets[0].indices()[0]
    frame.load_image(chooser_wrapper(imagesets[0],idx))
    frame.Show()
    app.MainLoop()
