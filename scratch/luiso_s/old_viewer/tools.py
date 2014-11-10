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

def show_reflection(reflection, orient = "landscape"):
  ''' Display a single reflection. '''
  import wx
  print "Show Reflection"

  class RefViewApp(wx.App):
    def OnInit(self):
      from dials.viewer.shoebox_view_frame import ShoeboxView
      self.frame = ShoeboxView(None, orient = orient, refl = reflection
                               , title = "Shoebox Viewer")
      self.SetTopWindow(self.frame)
      self.frame.Show()
      return True

  app = RefViewApp(redirect=False)
  app.MainLoop()
