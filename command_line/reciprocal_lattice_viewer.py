# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import absolute_import, division, print_function
from dials.util import wx_viewer
import copy
import wx
import wxtbx.utils
from gltbx.gl import *
import gltbx
from scitbx.math import minimum_covering_sphere
from dials.array_family import flex
from scitbx import matrix
import libtbx.phil
import dials.util.banner

help_message = '''

Visualise the strong spots from spotfinding in reciprocal space.

Examples::

  dials.reciprocal_lattice_viewer datablock.json strong.pickle

  dials.reciprocal_lattice_viewer experiments.json indexed.pickle

'''

phil_scope= libtbx.phil.parse("""
  reverse_phi = False
    .type = bool
    .optional = True
  beam_centre = None
    .type = floats(size=2)
    .help = "Fast, slow beam centre coordinates (mm)."
  show_rotation_axis = False
    .type = bool
  show_beam_vector = False
    .type = bool
  show_reciprocal_cell = False
    .type = bool
  d_min = None
    .type = float(value_min=0.0)
  z_min = None
    .type = float
  z_max = None
    .type = float
  n_min = 0
    .type = int
    .help = "Minimum size of spot in pixels"
  n_max = 0
    .type = int
    .help = "Maximum size of spot in pixels"
  partiality_min = None
    .type = float
  partiality_max = None
    .type = float
  display = *all unindexed indexed integrated
    .type = choice
  outlier_display = outliers inliers
    .type = choice
  marker_size = 5
    .type = int(value_min=1)
  experiment_ids = None
    .type = ints(value_min=-1)
  imageset_ids = None
    .type = ints(value_min=0)
  autospin = False
    .type = bool
  model_view_matrix = None
    .type = floats(size=16)
  black_background = True
    .type = bool
    .help = "Switch between black or white background"
""")

# WX4 compatibility
def _rewrite_event(unbound):
  """Decorator to intercept the event and add missing instance methods"""
  def _wrapp(self, event):
    event.GetPositionTuple = event.GetPosition
    return unbound(self, event)
  return _wrapp

# HACK: Monkeypatch wxtbx so that we don't use old interfaces
from wxtbx import segmentedctrl
segmentedctrl.SegmentedControl.HitTest = _rewrite_event(segmentedctrl.SegmentedControl.HitTest)
segmentedctrl.SegmentedControl.OnMotion = _rewrite_event(segmentedctrl.SegmentedControl.OnMotion)

def settings():
  return phil_scope.fetch().extract()

class render_3d(object):

  def __init__(self):
    self.reflections = None
    self.goniometer_orig = None

  def load_models(self, imagesets, reflections, crystals):
    self.imagesets = imagesets
    self.reflections_input = reflections
    if self.imagesets[0].get_goniometer() is not None:
      self.viewer.set_rotation_axis(
        self.imagesets[0].get_goniometer().get_rotation_axis())
    self.viewer.set_beam_vector(self.imagesets[0].get_beam().get_s0())
    detector = self.imagesets[0].get_detector()
    beam = self.imagesets[0].get_beam()
    if self.settings.beam_centre is None:
      try:
        panel_id, self.settings.beam_centre \
          = detector.get_ray_intersection(beam.get_s0())
      except RuntimeError:
        pass
    else:
      self.set_beam_centre(self.settings.beam_centre)
    if crystals is not None:
      vecs = [matrix.sqr(c.get_A()).transpose().as_list_of_lists() for c in crystals]
      self.viewer.set_reciprocal_lattice_vectors(vecs)
    self.map_points_to_reciprocal_space()
    self.set_points()

  def map_points_to_reciprocal_space(self):

    from dials.algorithms.indexing import indexer
    import copy

    reflections = flex.reflection_table()
    for i, imageset in enumerate(self.imagesets):
      if 'imageset_id' in self.reflections_input:
        sel = (self.reflections_input['imageset_id'] == i)
      else:
        sel = (self.reflections_input['id'] == i)
      if isinstance(self.reflections_input['id'], flex.size_t):
        self.reflections_input['id'] = self.reflections_input['id'].as_int()

      # 155 handle data from predictions *only* if that is what we have
      if 'xyzobs.px.value' in self.reflections_input:
        refl = indexer.indexer_base.map_spots_pixel_to_mm_rad(
          self.reflections_input.select(sel),
          imageset.get_detector(), imageset.get_scan())

        goniometer = copy.deepcopy(imageset.get_goniometer())
        if self.settings.reverse_phi:
          goniometer.set_rotation_axis(
            [-i for i in goniometer.get_rotation_axis()])
        indexer.indexer_base.map_centroids_to_reciprocal_space(
          refl, imageset.get_detector(), imageset.get_beam(),
          goniometer)

      else:
        # work on xyzcal.mm
        refl = self.reflections_input.select(sel)

        goniometer = copy.deepcopy(imageset.get_goniometer())
        if self.settings.reverse_phi:
          goniometer.set_rotation_axis(
            [-i for i in goniometer.get_rotation_axis()])

        indexer.indexer_base.map_centroids_to_reciprocal_space(
          refl, imageset.get_detector(), imageset.get_beam(),
          goniometer, calculated=True)

      reflections.extend(refl)
      self.reflections = reflections

  def set_points(self):
    reflections = self.reflections

    if 'miller_index' in reflections:
      if 'flags' not in reflections:
        reflections.set_flags(
          reflections['miller_index'] != (0,0,0), reflections.flags.indexed)
        reflections['id'].set_selected(
          reflections['miller_index'] == (0,0,0), -1)
        reflections.set_flags(
          flex.bool(len(reflections), True), reflections.flags.strong)
        reflections.set_flags(
          flex.bool(len(reflections), False), reflections.flags.integrated)
        reflections.set_flags(
          flex.bool(len(reflections), False), reflections.flags.centroid_outlier)

      outlier_sel = reflections.get_flags(reflections.flags.centroid_outlier)

      if self.settings.outlier_display == 'outliers':
        reflections = reflections.select(outlier_sel)
      if self.settings.outlier_display == 'inliers':
        reflections = reflections.select(~outlier_sel)

      indexed_sel = reflections.get_flags(reflections.flags.indexed)
      strong_sel = reflections.get_flags(reflections.flags.strong)
      integrated_sel = reflections.get_flags(reflections.flags.integrated)

      if self.settings.display == 'indexed':
        reflections = reflections.select(indexed_sel)
      elif self.settings.display == 'unindexed':
        reflections = reflections.select(strong_sel & ~indexed_sel)
      elif self.settings.display == 'integrated':
        reflections = reflections.select(integrated_sel)

      if self.settings.experiment_ids:
        sel = flex.bool(len(reflections), False)
        for i in self.settings.experiment_ids:
          sel.set_selected(reflections['id'] == i, True)
        reflections = reflections.select(sel)

    if self.settings.imageset_ids and 'imageset_id' in reflections:
      sel = flex.bool(len(reflections), False)
      for i in self.settings.imageset_ids:
        sel.set_selected(reflections['imageset_id'] == i, True)
      reflections = reflections.select(sel)

    d_spacings = 1/reflections['rlp'].norms()

    # 155 handle data from predictions *only* if that is what we have
    if 'xyzobs.px.value' in self.reflections_input:
      use_column = 'xyzobs.px.value'
    else:
      use_column = 'xyzcal.px'

    if self.settings.d_min is not None:
      reflections = reflections.select(d_spacings > self.settings.d_min)
    else:
      self.settings.d_min = flex.min(d_spacings)
    if self.settings.z_min is not None:
      z = reflections[use_column].parts()[2]
      reflections = reflections.select(z >= self.settings.z_min)
    else:
      z = reflections[use_column].parts()[2]
      self.settings.z_min = flex.min(z)
    if self.settings.z_max is not None:
      z = reflections[use_column].parts()[2]
      reflections = reflections.select(z <= self.settings.z_max)
    else:
      z = reflections[use_column].parts()[2]
      self.settings.z_max = flex.max(z)

    # N.B. 0 sensible NULL value here as cannot have 0 signal and still be spot

    if 'n_signal' in reflections:
      if self.settings.n_min != 0:
        _ns = reflections['n_signal']
        reflections = reflections.select(_ns >= self.settings.n_min)
      else:
        _ns = reflections['n_signal']
        self.settings.n_min = int(flex.min(_ns))

      if self.settings.n_max != 0:
        _ns = reflections['n_signal']
        reflections = reflections.select(_ns <= self.settings.n_max)
      else:
        _ns = reflections['n_signal']
        self.settings.n_max = int(flex.max(_ns))

    if 'partiality' in reflections:
      p = reflections['partiality']
      if self.settings.partiality_min is not None:
        reflections = reflections.select(p >= self.settings.partiality_min)
      else:
        self.settings.partiality_min = flex.min(p)
      if self.settings.partiality_max is not None:
        reflections = reflections.select(p <= self.settings.partiality_max)
      else:
        self.settings.partiality_max = flex.max(p)
    points = reflections['rlp'] * 100
    self.viewer.set_points(points)
    colors = flex.vec3_double(len(points), (1,1,1))

    if len(points):
      # suggested colorblind color palette
      # sorry if you have > 8 lattices!
      palette = flex.vec3_double((
        (255,255,255), (230,159,0), (86,180,233), (0,158,115),
        (240,228,66), (0,114,178), (213,94,0), (204,121,167)))
      palette *= (1/255)
      if not self.viewer.settings.black_background:
        bkg = flex.vec3_double()
        for j in range(len(palette)):
          bkg.append((1,1,1))
        palette = bkg - palette
      self.viewer.set_palette(palette)
      n = palette.size() - 1
      if self.reflections.get_flags(self.reflections.flags.indexed).count(True) == 0:
        if 'imageset_id' in reflections:
          imageset_id = reflections['imageset_id']
        else:
          imageset_id = reflections['id']
        for i in range(0, flex.max(imageset_id)+1):
          colors.set_selected(imageset_id == i, palette[(i%n)+1])
      else:
        colors.set_selected(reflections['id'] == -1, palette[0])
        for i in range(0, flex.max(reflections['id'])+1):
          colors.set_selected(reflections['id'] == i, palette[(i%n)+1])
    self.viewer.set_colors(colors)


class ReciprocalLatticeViewer(wx.Frame, render_3d):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    render_3d.__init__(self)
    self.parent = self.GetParent()
    self.statusbar = self.CreateStatusBar()
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)

    app = wx.GetApp()
    if getattr(app, "settings", None) is not None:
      # XXX copying the initial settings avoids awkward interactions when
      # multiple viewer windows are opened
      self.settings = copy.deepcopy(app.settings)
    else :
      self.settings = settings()

    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND|wx.ALL)
    self.SetSizerAndFit(self.sizer)
    self.SetMinSize(self.settings_panel.GetSize())
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.Bind(wx.EVT_ACTIVATE, self.OnActive)
    self.viewer.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
    self.viewer.SetFocus()

  def load_models(self, imagesets, reflections, crystals):
    render_3d.load_models(self, imagesets, reflections, crystals)
    if self.settings.beam_centre is not None:
      self.settings_panel.beam_fast_ctrl.SetValue(self.settings.beam_centre[0])
      self.settings_panel.beam_slow_ctrl.SetValue(self.settings.beam_centre[1])
    self.settings_panel.marker_size_ctrl.SetValue(self.settings.marker_size)
    self.settings_panel.add_experiments_buttons()
    self.settings_panel.add_imagesets_buttons()

  def OnActive (self, event) :
    if self.IsShown() and type(self.viewer).__name__ != "_wxPyDeadObject":
      self.viewer.Refresh()

  def OnClose (self, event) :
    self.Unbind(wx.EVT_ACTIVATE)
    self.Destroy()
    event.Skip()

  def OnDestroy (self, event) :
    if self.parent is not None:
      self.parent.viewer = None
    event.Skip()

  def OnKeyDown(self, event):
    key = event.GetUnicodeKey()
    if key == wx.WXK_NONE:
      key = event.GetKeyCode()
    dxs = {wx.WXK_LEFT:-1,
           wx.WXK_RIGHT:+1,
           wx.WXK_UP:0,
           wx.WXK_DOWN:0}
    dys = {wx.WXK_LEFT:0,
           wx.WXK_RIGHT:0,
           wx.WXK_UP:+1,
           wx.WXK_DOWN:-1}

    if key in dxs:
      dx = dxs[key]
      dy = dys[key]
      if event.ShiftDown():
        scale = 0.1
      else:
        scale = 1.0
      self.do_Step(dx, dy, scale)

  def do_Step(self, dx, dy, scale):
    v = self.viewer
    rc = v.rotation_center
    glMatrixMode(GL_MODELVIEW)
    gltbx.util.rotate_object_about_eye_x_and_y(
      scale, rc[0], rc[1], rc[2],
      dx, dy, 0, 0)
    v.OnRedraw()

  def create_viewer_panel (self) :
    if self.settings.black_background:
      background_rgb = (0,0,0)
    else:
      background_rgb = (255,255,255)
    self.viewer = RLVWindow(
      settings=self.settings, parent=self, size=(800,600),
      background_rgb=background_rgb)

  def create_settings_panel (self) :
    self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

  def set_points(self):
    render_3d.set_points(self)
    self.settings_panel.d_min_ctrl.SetValue(self.settings.d_min)
    self.settings_panel.z_min_ctrl.SetValue(self.settings.z_min)
    self.settings_panel.z_max_ctrl.SetValue(self.settings.z_max)
    self.settings_panel.n_min_ctrl.SetValue(self.settings.n_min)
    self.settings_panel.n_max_ctrl.SetValue(self.settings.n_max)
    if self.settings.partiality_min is not None:
      self.settings_panel.partiality_min_ctrl.SetValue(self.settings.partiality_min)
    if self.settings.partiality_max is not None:
      self.settings_panel.partiality_max_ctrl.SetValue(self.settings.partiality_max)

  def update_settings(self, *args, **kwds):

    imageset = self.imagesets[0]
    detector = imageset.get_detector()
    beam = imageset.get_beam()

    try:
      panel_id, beam_centre = detector.get_ray_intersection(beam.get_s0())
    except RuntimeError:
      # beam centre calculation fails if the beam falls between panels
      pass
    else:
      if self.settings.beam_centre != beam_centre:
        self.set_beam_centre(beam_centre)

    self.map_points_to_reciprocal_space()
    self.set_points()
    self.viewer.update_settings(*args, **kwds)

  def set_beam_centre(self, beam_centre):

    imageset = self.imagesets[0]
    detector = imageset.get_detector()
    beam = imageset.get_beam()

    try:
      panel_id, beam_centre = detector.get_ray_intersection(beam.get_s0())
    except RuntimeError:
      # beam centre calculation fails if the beam falls between panels
      pass
    else:
      # code copied from: dials.command_line.dials_import.PixelBeamCenterUpdater

      beam_f, beam_s = self.settings.beam_centre

      try:
        p = detector[panel_id]
      except RuntimeError:
        raise Sorry('Detector does not have panel index {0}'.format(panel_id))

      beam_f, beam_s = p.millimeter_to_pixel((beam_f, beam_s))
      beam.set_unit_s0(p.get_pixel_lab_coord((beam_f, beam_s)))

  def update_statusbar (self) :
    model_view_matrix = gltbx.util.get_gl_modelview_matrix()
    txt = 'Model view matrix: ' + '[' + ', '.join(
      '%.4f' %m for m in model_view_matrix) + ']'
    self.statusbar.SetStatusText(txt)

class settings_window (wxtbx.utils.SettingsPanel) :
  def __init__ (self, *args, **kwds) :
    wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_CHAR, self.OnChar)

  def OnChar (self, event) :
    self.GetParent().viewer.OnChar(event)

  def add_controls (self) :
    # d_min control
    from wx.lib.agw import floatspin
    self.d_min_ctrl = floatspin.FloatSpin(parent=self, increment=0.05, digits=2)
    self.d_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.d_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"High resolution:")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.d_min_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.d_min_ctrl)

    self.z_min_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=0)
    self.z_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.z_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Min Z")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.z_min_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.z_min_ctrl)

    self.z_max_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=0)
    self.z_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.z_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Max Z")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.z_max_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.z_max_ctrl)

    # Control for spot size (utility depends on n_signal column in reflection
    # file - will be ignored if not in file

    self.n_min_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=0)
    self.n_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.n_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Min Pixels")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.n_min_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.n_min_ctrl)

    self.n_max_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=0)
    self.n_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.n_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Max Pixels")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.n_max_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.n_max_ctrl)

    # end new control


    self.partiality_min_ctrl = floatspin.FloatSpin(
      parent=self, increment=0.01, digits=3, min_val=0, max_val=1)
    self.partiality_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.partiality_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Min partiality")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.partiality_min_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.partiality_min_ctrl)

    self.partiality_max_ctrl = floatspin.FloatSpin(
      parent=self, increment=0.01, digits=3, min_val=0, max_val=1)
    self.partiality_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.partiality_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Max partiality")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.partiality_max_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.partiality_max_ctrl)

    ctrls = self.create_controls(
      setting="show_rotation_axis",
      label="Show rotation axis")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_beam_vector",
      label="Show beam vector")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_reciprocal_cell",
      label="Show reciprocal cell")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    self.reverse_phi_ctrl = self.create_controls(
      setting="reverse_phi",
      label="Invert rotation axis")[0]
    self.panel_sizer.Add(self.reverse_phi_ctrl, 0, wx.ALL, 5)
    self.Bind(wx.EVT_CHECKBOX, self.OnChangeSettings,
              self.reverse_phi_ctrl)

    self.beam_fast_ctrl = floatspin.FloatSpin(parent=self, increment=0.01,
                                              digits=2)
    self.beam_fast_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.beam_fast_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Beam centre (mm)")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.beam_fast_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings,
              self.beam_fast_ctrl)

    self.beam_slow_ctrl = floatspin.FloatSpin(parent=self, increment=0.01,
                                              digits=2)
    self.beam_slow_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.beam_slow_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box.Add(self.beam_slow_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.beam_slow_ctrl)

    self.marker_size_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=0,
                                          min_val=1)
    self.marker_size_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.marker_size_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Marker size:")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.marker_size_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.marker_size_ctrl)

    from wxtbx.segmentedctrl import SegmentedRadioControl, SEGBTN_HORIZONTAL
    self.btn = SegmentedRadioControl(self, style=SEGBTN_HORIZONTAL)
    self.btn.AddSegment("all")
    self.btn.AddSegment("indexed")
    self.btn.AddSegment("unindexed")
    self.btn.AddSegment("integrated")
    self.btn.SetSelection(
      ["all", "indexed", "unindexed", "integrated"].index(self.settings.display))
    self.Bind(wx.EVT_RADIOBUTTON, self.OnChangeSettings, self.btn)
    self.GetSizer().Add(self.btn, 0, wx.ALL, 5)

    self.outlier_btn = SegmentedRadioControl(self, style=SEGBTN_HORIZONTAL)
    self.outlier_btn.AddSegment("all")
    self.outlier_btn.AddSegment("inliers")
    self.outlier_btn.AddSegment("outliers")
    self.outlier_btn.SetSelection(
      [None, "inliers", "outliers"].index(self.settings.outlier_display))
    self.Bind(wx.EVT_RADIOBUTTON, self.OnChangeSettings, self.outlier_btn)
    self.GetSizer().Add(self.outlier_btn, 0, wx.ALL, 5)

  def add_value_widgets (self, sizer) :
    sizer.Add(wx.StaticText(self.panel, -1, "Value:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.value_info = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    sizer.Add(self.value_info, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def add_experiments_buttons(self):
    n = flex.max(self.parent.reflections_input['id'])
    if n <= 0:
      self.expt_btn = None
      return

    box = wx.BoxSizer(wx.VERTICAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Experiment ids:")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    from wxtbx.segmentedctrl import SegmentedToggleControl, SEGBTN_HORIZONTAL
    self.expt_btn = SegmentedToggleControl(self, style=SEGBTN_HORIZONTAL)
    for i in range(-1, n+1):
      self.expt_btn.AddSegment(str(i))
      if (self.settings.experiment_ids is not None and
          i in self.settings.experiment_ids):
        self.expt_btn.SetValue(i+1, True)

    self.expt_btn.Realize()
    self.Bind(wx.EVT_TOGGLEBUTTON, self.OnChangeSettings, self.expt_btn)
    box.Add(self.expt_btn, 0, wx.ALL, 5)

  def add_imagesets_buttons(self):
    if 'imageset_id' not in self.parent.reflections_input:
      self.imgset_btn = None
      return
    n = flex.max(self.parent.reflections_input['imageset_id'])
    if n <= 0:
      self.imgset_btn = None
      return

    box = wx.BoxSizer(wx.VERTICAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"Imageset ids:")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    from wxtbx.segmentedctrl import SegmentedToggleControl, SEGBTN_HORIZONTAL
    self.imgset_btn = SegmentedToggleControl(self, style=SEGBTN_HORIZONTAL)
    for i in range(n+1):
      self.imgset_btn.AddSegment(str(i))
      if (self.settings.imageset_ids is not None and
          i in self.settings.imageset_ids):
        self.imgset_btn.SetValue(i+1, True)

    self.imgset_btn.Realize()
    self.Bind(wx.EVT_TOGGLEBUTTON, self.OnChangeSettings, self.imgset_btn)
    box.Add(self.imgset_btn, 0, wx.ALL, 5)

  def OnChangeSettings(self, event):
    self.settings.d_min = self.d_min_ctrl.GetValue()
    self.settings.z_min = self.z_min_ctrl.GetValue()
    self.settings.z_max = self.z_max_ctrl.GetValue()
    self.settings.n_min = int(self.n_min_ctrl.GetValue())
    self.settings.n_max = int(self.n_max_ctrl.GetValue())
    self.settings.partiality_min = self.partiality_min_ctrl.GetValue()
    self.settings.partiality_max = self.partiality_max_ctrl.GetValue()
    self.settings.beam_centre = (
      self.beam_fast_ctrl.GetValue(), self.beam_slow_ctrl.GetValue())
    self.settings.reverse_phi = self.reverse_phi_ctrl.GetValue()
    self.settings.marker_size = self.marker_size_ctrl.GetValue()
    for i, display in enumerate(("all", "indexed", "unindexed", "integrated")):
      if self.btn.values[i]:
        self.settings.display = display
        break

    for i, display in enumerate(("all", "inliers", "outliers")):
      if self.outlier_btn.values[i]:
        self.settings.outlier_display = display
        break

    if self.expt_btn is not None:
      expt_ids = []
      for i in range(len(self.expt_btn.segments)):
        if self.expt_btn.GetValue(i):
          expt_ids.append(i-1)
      self.settings.experiment_ids = expt_ids

    if self.imgset_btn is not None:
      imgset_ids = []
      for i in range(len(self.imgset_btn.segments)):
        if self.imgset_btn.GetValue(i):
          imgset_ids.append(i)
      self.settings.imageset_ids = imgset_ids

    self.parent.update_settings()


class RLVWindow(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, settings, *args, **kwds):
    super(RLVWindow, self).__init__(*args, **kwds)
    self.settings = settings
    self.points = flex.vec3_double()
    self.colors = None
    self.palette = None
    self.rotation_axis = None
    self.beam_vector = None
    self.recip_latt_vectors = None
    self.flag_show_minimum_covering_sphere = False
    self.minimum_covering_sphere = None
    self.field_of_view_y = 0.001
    if self.settings.autospin:
      self.autospin_allowed = True
      self.yspin = 1
      self.xspin = 1
      self.autospin = True

  def set_points(self, points):
    self.points = points
    self.points_display_list = None
    if self.minimum_covering_sphere is None:
      self.update_minimum_covering_sphere()

  def set_colors(self, colors):
    assert len(colors) == len(self.points)
    self.colors = colors

  def set_palette(self, palette):
    self.palette = palette

  def draw_points(self):
    if self.points_display_list is None:
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      glLineWidth(1)
      if self.colors is None:
        self.colors = flex.vec3_double(len(self.points), (1,1,1))
      for point, color in zip(self.points, self.colors):
        self.draw_cross_at(point, color=color)
      self.points_display_list.end()
    self.points_display_list.call()

  def set_rotation_axis(self, axis):
    self.rotation_axis = axis

  def set_beam_vector(self, beam):
    self.beam_vector = beam

  def set_reciprocal_lattice_vectors(self, vectors_per_crystal):
    self.recip_latt_vectors = []
    # the points are scaled by 100 so must do that here too
    for vec_set in vectors_per_crystal:
      self.recip_latt_vectors.append([
        100. * matrix.col(vec) for vec in vec_set])

  #--- user input and settings
  def update_settings (self) :
    self.points_display_list = None
    self.Refresh()

  def update_minimum_covering_sphere(self):
    n_points = min(1000, self.points.size())
    isel = flex.random_permutation(self.points.size())[:n_points]
    self.minimum_covering_sphere = minimum_covering_sphere(
      self.points.select(isel))

  def draw_cross_at(self, xyz, color=(1,1,1), f=None):
    (x,y,z) = xyz
    if f is None:
      f = 0.01 * self.settings.marker_size
    wx_viewer.show_points_and_lines_mixin.draw_cross_at(
      self, (x,y,z), color=color, f=f)

  def DrawGL(self):
    wx_viewer.show_points_and_lines_mixin.DrawGL(self)
    if self.rotation_axis is not None and self.settings.show_rotation_axis:
      self.draw_axis(self.rotation_axis, "phi")
    if self.beam_vector is not None and self.settings.show_beam_vector:
      self.draw_axis(self.beam_vector, "beam")
    if self.recip_latt_vectors is not None and self.settings.show_reciprocal_cell:
      for i, axes in enumerate(self.recip_latt_vectors):
        if self.settings.experiment_ids:
          if i not in self.settings.experiment_ids: continue
        j = (i + 1) % self.palette.size()
        color = self.palette[j]
        self.draw_cell(axes, color)
    self.GetParent().update_statusbar()

  def draw_axis(self, axis, label):
    if self.minimum_covering_sphere is None:
      self.update_minimum_covering_sphere()
    s = self.minimum_covering_sphere
    scale = max(max(s.box_max()), abs(min(s.box_min())))
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    glDisable(GL_LIGHTING)
    glColor3f(1.0, 1.0, 1.0)
    glLineWidth(1.0)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(axis[0]*scale, axis[1]*scale, axis[2]*scale)
    glEnd()
    glRasterPos3f(0.5+axis[0]*scale, 0.2+axis[1]*scale, 0.2+axis[2]*scale)
    gltbx.fonts.ucs_bitmap_8x13.render_string(label)
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-axis[0]*scale, -axis[1]*scale, -axis[2]*scale)
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  def draw_cell(self, axes, color):
    astar, bstar, cstar = axes[0], axes[1], axes[2]
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    glDisable(GL_LIGHTING)
    glColor3f(*color)
    glLineWidth(2.0)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(*astar.elems)
    glVertex3f(0.,0.,0.)
    glVertex3f(*bstar.elems)
    glVertex3f(0.,0.,0.)
    glVertex3f(*cstar.elems)
    glEnd()
    glRasterPos3f(*(1.01 * astar).elems)
    gltbx.fonts.ucs_bitmap_8x13.render_string("a*")
    glRasterPos3f(*(1.01 * bstar).elems)
    gltbx.fonts.ucs_bitmap_8x13.render_string("b*")
    glRasterPos3f(*(1.01 * cstar).elems)
    gltbx.fonts.ucs_bitmap_8x13.render_string("c*")
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    farpoint = astar + bstar + cstar
    # a* face
    glBegin(GL_LINE_LOOP)
    glVertex3f(*farpoint.elems)
    glVertex3f(*(farpoint - bstar).elems)
    glVertex3f(*(farpoint - bstar - cstar).elems)
    glVertex3f(*(farpoint - cstar).elems)
    glEnd()
    # b* face
    glBegin(GL_LINE_LOOP)
    glVertex3f(*farpoint.elems)
    glVertex3f(*(farpoint - astar).elems)
    glVertex3f(*(farpoint - astar - cstar).elems)
    glVertex3f(*(farpoint - cstar).elems)
    glEnd()
    # c* face
    glBegin(GL_LINE_LOOP)
    glVertex3f(*farpoint.elems)
    glVertex3f(*(farpoint - bstar).elems)
    glVertex3f(*(farpoint - bstar - astar).elems)
    glVertex3f(*(farpoint - astar).elems)
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  def rotate_view(self, x1, y1, x2, y2, shift_down=False, scale=0.1):
    super(RLVWindow, self).rotate_view(
      x1, y1, x2, y2, shift_down=shift_down, scale=scale)

  def OnLeftUp(self,event):
    self.was_dragged = True
    super(RLVWindow, self).OnLeftUp(event)

  def initialize_modelview(self, eye_vector=None, angle=None):
    super(RLVWindow, self).initialize_modelview(eye_vector=eye_vector, angle=angle)
    self.rotation_center = (0,0,0)
    self.move_to_center_of_viewport(self.rotation_center)
    if self.settings.model_view_matrix is not None:
      glLoadMatrixd(self.settings.model_view_matrix)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json reflections.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if (len(datablocks) == 0 and len(experiments) == 0) or len(reflections) == 0:
    parser.print_help()
    exit(0)

  if len(datablocks) == 0 and len(experiments) > 0:
    imagesets = experiments.imagesets()
    crystals = experiments.crystals()
  else:
    imagesets = []
    crystals = None
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  if len(reflections) > 1:
    assert len(reflections) == len(imagesets)
    for i in range(len(reflections)):
      reflections[i]['imageset_id'] = flex.int(len(reflections[i]), i)
      if i > 0:
        reflections[0].extend(reflections[i])
  elif 'imageset_id' not in reflections[0]:
    reflections[0]['imageset_id'] = reflections[0]['id']
    reflections[0]['id'] = flex.int(reflections[0].size(), -1)

  reflections = reflections[0]

  import wxtbx.app
  a = wxtbx.app.CCTBXApp(0)
  a.settings = params
  f = ReciprocalLatticeViewer(
    None, -1, "Reflection data viewer", size=(1024,768))
  f.load_models(imagesets, reflections, crystals)
  f.Show()
  a.SetTopWindow(f)
  #a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
