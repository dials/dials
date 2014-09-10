#
#  DIALS Image Utilities
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luiso and James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

from __future__ import division

def show_reflection(reflection, **kwargs):
  ''' Display a single reflection. '''
  import wx
  print "Show Reflection"

  class RefViewApp(wx.App):
    def OnInit(self):
      from dials.viewer.image_frame import ImageFrame
      self.frame = ImageFrame(None, refl = reflection, title = "Reflection")
      self.SetTopWindow(self.frame)
      self.frame.Show()
      return True

  app = RefViewApp(redirect=False)
  app.MainLoop()
