from __future__ import division
class spot_wrapper:
  def __init__(self,working_phil):
    self.working_phil = working_phil

  def display(self, sweep_filenames, reflections):
    import wx
    from dials.util.spotfinder_frame import SpotFrame

    app   = wx.App(0)
    frame = SpotFrame(None, -1, "X-ray image display", size=(800,720),
      pos=(100,100),
      sweep_filenames=sweep_filenames,
      reflections=reflections)
    frame.SetSize((1024,780))

    file_paths = sweep_filenames
    if len(file_paths) > 0:
      from dxtbx.imageset import ImageSetFactory
      import time
      t = time.time()
      print "Starting loading files"
      imgsets = ImageSetFactory.new(file_paths, ignore_unknown=True)
      print "Finished loading images (time taken = %s)" %(time.time() - t)

    from rstbx.slip_viewer.frame import chooser_wrapper

    i = 0
    for imgset in imgsets:
      for idx in imgset.indices():
        i += 1
        print i
        frame.add_file_name_or_data(chooser_wrapper(imgset, idx))
    idx = imgsets[0].indices()[0]
    frame.load_image(chooser_wrapper(imgsets[0],idx))
    frame.Show()
    app.MainLoop()
