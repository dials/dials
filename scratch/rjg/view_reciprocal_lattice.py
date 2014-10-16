from __future__ import division
from gltbx import wx_viewer
import copy
import wx
import wxtbx.utils
from gltbx.gl import *
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
import libtbx.phil

master_phil = libtbx.phil.parse("""
  data = None
    .type = path
    .optional = False
  reverse_phi = False
    .type = bool
    .optional = True
""")

def settings () :
  return master_phil.fetch().extract()

class ReciprocalLatticeViewer(wx.Frame):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    self.parent = self.GetParent()
    self.statusbar = self.CreateStatusBar()
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.detector = None
    self.beam = None
    self.scan = None
    self.goniometer = None
    self.reflections = None

    app = wx.GetApp()
    if (getattr(app, "settings", None) is not None) :
      # XXX copying the initial settings avoids awkward interactions when
      # multiple viewer windows are opened
      self.settings = copy.deepcopy(app.settings)
    else :
      self.settings = settings()

    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND|wx.ALL)
    #self.SetupToolbar()
    #self.SetupMenus()
    #self.add_view_specific_functions()
    #self.SetMenuBar(self.menubar)
    #self.toolbar.Realize()
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.Bind(wx.EVT_ACTIVATE, self.OnActive)
    self.viewer.SetFocus()

  def OnActive (self, event) :
    if (self.IsShown()) and (type(self.viewer).__name__ != "_wxPyDeadObject") :
      self.viewer.Refresh()

  def OnClose (self, event) :
    self.Unbind(wx.EVT_ACTIVATE)
    self.Destroy()
    event.Skip()

  def OnDestroy (self, event) :
    if (self.parent is not None) :
      self.parent.viewer = None
    event.Skip()

  def create_viewer_panel (self) :
    self.viewer = MyGLWindow(parent=self, size=(800,600),
      style=wx.glcanvas.WX_GL_DOUBLEBUFFER,
      #orthographic=True
      )

  def create_settings_panel (self) :
    self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

  def load_models(self, imageset, reflections):
    self.detector = imageset.get_detector()
    self.beam = imageset.get_beam()
    self.scan = imageset.get_scan()
    self.goniometer = imageset.get_goniometer()
    self.reflections = reflections
    self.parameterise_models()
    self.map_points_to_reciprocal_space()

  def parameterise_models(self):
    from dials.algorithms.refinement.parameterisation import \
      DetectorParameterisationSinglePanel, DetectorParameterisationHierarchical, \
      DetectorParameterisationMultiPanel
    if len(self.detector) > 1:
      try:
        h = self.detector.hierarchy()
        self.dp = DetectorParameterisationHierarchical(self.detector,
            experiment_ids=[0], level=0)
      except AttributeError:
        self.dp = DetectorParameterisationMultiPanel(self.detector, self.beam,
                                                    experiment_ids=[0])
    else:
      self.dp = DetectorParameterisationSinglePanel(self.detector,
                                                    experiment_ids=[0])

    from dials.algorithms.refinement.parameterisation import BeamParameterisation
    self.bp = BeamParameterisation(self.beam, self.goniometer, experiment_ids=[0])

    # keep parameter names
    self.detector_parameter_names = self.dp.get_param_names()
    self.beam_parameter_names = self.bp.get_param_names()

    # keep parameter values (modify these and update_models to make changes)
    self.detector_parameters = self.dp.get_param_vals()
    self.beam_parameters = self.bp.get_param_vals()

  def map_points_to_reciprocal_space(self):
    goniometer = copy.deepcopy(self.goniometer)
    if self.settings.reverse_phi:
      goniometer.set_rotation_axis([-i for i in goniometer.get_rotation_axis()])
    from dials.algorithms.indexing import indexer

    # set parameters to latest values
    self.dp.set_param_vals(self.detector_parameters)
    self.bp.set_param_vals(self.beam_parameters)

    reflections = indexer.indexer_base.map_spots_pixel_to_mm_rad(
      self.reflections, self.detector, self.scan)
    indexer.indexer_base.map_centroids_to_reciprocal_space(
      reflections, self.detector, self.beam, goniometer)
    points = reflections['rlp'] * 100
    self.viewer.set_points(points)

  def update_settings(self, *args, **kwds):
    self.map_points_to_reciprocal_space()
    self.viewer.update_settings(*args, **kwds)


class settings_window (wxtbx.utils.SettingsPanel) :
  def __init__ (self, *args, **kwds) :
    wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_CHAR, self.OnChar)

  def OnChar (self, event) :
    self.GetParent().viewer.OnChar(event)

  def add_controls (self) :
    ctrls = self.create_controls(
      setting="reverse_phi",
      label="Reverse phi direction")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)

  def add_value_widgets (self, sizer) :
    sizer.Add(wx.StaticText(self.panel, -1, "Value:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.value_info = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    sizer.Add(self.value_info, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def OnChangeColor (self, event) :
    self.settings.color_scheme = str(self.color_ctrl.GetStringSelection())
    self.parent.update_settings()

  def OnSetScale (self, event) :
    self.settings.scale = (self.scale_ctrl.GetValue() + 4) / 4
    self.parent.update_settings()


class MyGLWindow(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, *args, **kwds):
    super(MyGLWindow, self).__init__(*args, **kwds)
    self.points = flex.vec3_double()
    self.flag_show_minimum_covering_sphere = False
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)

  def set_points(self, points):
    self.points = points
    self.points_display_list = None
    #self.draw_points()
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)
    if not self.GL_uninitialised:
      self.fit_into_viewport()

  #--- user input and settings
  def update_settings (self) :
    self.points_display_list = None
    self.draw_points()
    self.minimum_covering_sphere = minimum_covering_sphere(self.points)
    if not self.GL_uninitialised:
      self.fit_into_viewport()
    self.Refresh()


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_reflections

  parser = OptionParser(
    phil=master_phil,
    read_datablocks=True,
    read_reflections=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)[0]
  assert len(datablocks) == 1
  imageset = datablocks[0].extract_imagesets()[0]

  import wxtbx.app
  a = wxtbx.app.CCTBXApp(0)
  a.settings = params
  f = ReciprocalLatticeViewer(
    None, -1, "Reflection data viewer", size=(1024,768))
  f.load_models(imageset, reflections)
  f.Show()
  a.SetTopWindow(f)
  #a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
