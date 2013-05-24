from __future__ import division
class spot_wrapper:
  def __init__(self,working_phil):
      self.working_phil = working_phil

  def display(self, sweep_filenames, reflections):
    import wx
    from dials.util.spotfinder_frame import SpotFrame
    from rstbx.viewer import display

    app   = wx.App(0)
    frame = SpotFrame(None, -1, "X-ray image display", size=(1200,1080),
      pos=(100,100),
      sweep_filenames=sweep_filenames,
      reflections=reflections)
    frame.SetSize((1024,780))
    for filename in sweep_filenames:
      frame.add_file_name_or_data(filename)
    path = sweep_filenames[0]
    frame.load_image(path)
    frame.path = path
    self.path = path
    frame.Show()
    app.MainLoop()
