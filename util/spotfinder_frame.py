from __future__ import division
import rstbx.viewer.display
import wx
from rstbx.slip_viewer.frame import XrayFrame
from rstbx.viewer.frame import SettingsFrame, SettingsPanel

class SpotFrame(XrayFrame) :
  def __init__ (self, *args, **kwds) :
    self.imagesets = kwds["imagesets"]
    self.reflections = kwds["reflections"]
    self.crystals = kwds["crystals"]
    del kwds["imagesets"]; del kwds["reflections"] #otherwise wx complains
    del kwds["crystals"] #otherwise wx complains
    super(SpotFrame, self).__init__(*args, **kwds)
    self.viewer.reflections = self.reflections
    self.viewer.frames = self.imagesets
    self.dials_spotfinder_layer = None
    self.shoebox_layer = None
    self.ctr_mass_layer = None
    self.max_pix_layer = None
    self.predictions_layer = None
    self.miller_indices_layer = None
    self.vector_layer = None
    self.vector_text_layer = None

    from libtbx.utils import time_log
    self.show_all_pix_timer = time_log("show_all_pix")
    self.show_shoebox_timer = time_log("show_shoebox")
    self.show_max_pix_timer = time_log("show_max_pix")
    self.show_ctr_mass_timer = time_log("show_ctr_mass")
    self.draw_all_pix_timer = time_log("draw_all_pix")
    self.draw_shoebox_timer = time_log("draw_shoebox")
    self.draw_max_pix_timer = time_log("draw_max_pix")
    self.draw_ctr_mass_timer = time_log("draw_ctr_mass_pix")

  #def __del__(self):
    #print self.show_all_pix_timer.legend
    #print self.show_all_pix_timer.report()
    #print self.show_shoebox_timer.report()
    #print self.show_max_pix_timer.report()
    #print self.show_ctr_mass_timer.report()
    #print self.draw_all_pix_timer.report()
    #print self.draw_shoebox_timer.report()
    #print self.draw_max_pix_timer.report()
    #print self.draw_ctr_mass_timer.report()

  def OnShowSettings (self, event) :
    if (self.settings_frame is None) :
      frame_rect = self.GetRect()
      display_rect = wx.GetClientDisplayRect()
      x_start = frame_rect[0] + frame_rect[2]
      if (x_start > (display_rect[2] - 400)) :
        x_start = display_rect[2] - 400
      y_start = frame_rect[1]
      self.settings_frame = SpotSettingsFrame(self, -1, "Settings",
        style=wx.CAPTION|wx.MINIMIZE_BOX, pos=(x_start, y_start))
    self.settings_frame.Show()

  def show_filters(self):
    from dials.algorithms.image.threshold import KabschDebug
    from dials.array_family import flex

    image = self.pyslip.tiles.raw_image
    raw_data = image.get_raw_data()
    detector = image.get_detector()

    trange = [p.get_trusted_range() for p in detector]
    mask = []
    for tr, im in zip(trange, [raw_data]):
      mask.append(im > int(tr[0]))

    nsigma_b = self.settings.nsigma_b
    nsigma_s = self.settings.nsigma_s
    min_count = 2
    size = self.settings.kernel_size
    debug = KabschDebug(image.get_raw_data().as_double(),
      mask[0], size, nsigma_b, nsigma_s, min_count)
    mean = debug.mean()
    variance = debug.variance()
    cv = debug.coefficient_of_variation()
    cv_mask = debug.cv_mask()
    value_mask = debug.value_mask()
    final_mask = debug.final_mask()
    if self.settings.show_mean_filter:
      self.pyslip.tiles.set_image_data(mean)
    elif self.settings.show_variance_filter:
      self.pyslip.tiles.set_image_data(variance)
    elif self.settings.show_dispersion:
      self.pyslip.tiles.set_image_data(cv)
    elif self.settings.show_sigma_b_filter:
      cv_mask = cv_mask.as_1d().as_double()
      cv_mask.reshape(mean.accessor())
      self.pyslip.tiles.set_image_data(cv_mask)
    elif self.settings.show_sigma_s_filter:
      value_mask = value_mask.as_1d().as_double()
      value_mask.reshape(mean.accessor())
      self.pyslip.tiles.set_image_data(value_mask)
    elif self.settings.show_threshold_map:
      final_mask = final_mask.as_1d().as_double()
      final_mask.reshape(mean.accessor())
      self.pyslip.tiles.set_image_data(final_mask)
    else:
      self.pyslip.tiles.set_image_data(raw_data)
    self.pyslip.ZoomToLevel(self.pyslip.tiles.zoom_level)
    self.update_statusbar() # XXX Not always working?
    self.Layout()

  def update_settings(self, layout=True):
    super(SpotFrame, self).update_settings(layout=layout)
    if self.settings.show_dials_spotfinder_spots:
      spotfinder_data = self.get_spotfinder_data()
      shoebox_data = spotfinder_data.shoebox_data
      all_pix_data = spotfinder_data.all_pix_data
      ctr_mass_data = spotfinder_data.ctr_mass_data
      max_pix_data = spotfinder_data.max_pix_data
      predictions_data = spotfinder_data.predictions_data
      miller_indices_data = spotfinder_data.miller_indices_data
      vector_data = spotfinder_data.vector_data
      vector_text_data = spotfinder_data.vector_text_data
      if self.dials_spotfinder_layer is not None:
        self.pyslip.DeleteLayer(self.dials_spotfinder_layer)
        self.dials_spotfinder_layer = None
      if self.shoebox_layer is not None:
        self.pyslip.DeleteLayer(self.shoebox_layer)
        self.shoebox_layer = None
      if self.ctr_mass_layer is not None:
        self.pyslip.DeleteLayer(self.ctr_mass_layer)
        self.ctr_mass_layer = None
      if self.max_pix_layer is not None:
        self.pyslip.DeleteLayer(self.max_pix_layer)
        self.max_pix_layer = None
      if self.predictions_layer is not None:
        self.pyslip.DeleteLayer(self.predictions_layer)
        self.predictions_layer = None
      if self.miller_indices_layer is not None:
        self.pyslip.DeleteLayer(self.miller_indices_layer)
        self.miller_indices_layer = None
      if self.vector_layer is not None:
        self.pyslip.DeleteLayer(self.vector_layer)
        self.vector_layer = None
      if self.vector_text_layer is not None:
        self.pyslip.DeleteLayer(self.vector_text_layer)
        self.vector_text_layer = None

      if self.settings.show_miller_indices and len(miller_indices_data):
        self.miller_indices_layer = self.pyslip.AddTextLayer(
          miller_indices_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<miller_indices_layer>')
      if self.settings.show_predictions and len(predictions_data):
        self.predictions_layer = self.pyslip.AddPointLayer(
          predictions_data, color="yellow", name="<predictions_layer>",
          radius=3,
          renderer = self.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      if self.settings.show_all_pix:
        self.draw_all_pix_timer.start()
        self.dials_spotfinder_layer = self.pyslip.AddPointLayer(
          all_pix_data, color="green", name="<all_pix_layer>",
          radius=2,
          renderer = self.pyslip.LightweightDrawPointLayer2,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
        self.draw_all_pix_timer.stop()
      if self.settings.show_shoebox:
        self.draw_shoebox_timer.start()
        self.shoebox_layer = self.pyslip.AddPolygonLayer(
          shoebox_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<shoebox_layer>')
        self.draw_shoebox_timer.stop()
      if self.settings.show_ctr_mass:
        self.draw_ctr_mass_timer.start()
        self.ctr_mass_layer = self.pyslip.AddPolygonLayer(
          ctr_mass_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<ctr_mass_layer>')
        self.draw_ctr_mass_timer.stop()
      if self.settings.show_max_pix:
        self.draw_max_pix_timer.start()
        self.max_pix_layer = self.pyslip.AddPointLayer(
          max_pix_data, color="pink", name="<max_pix_layer>",
          radius=2,
          renderer = self.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
        self.draw_max_pix_timer.stop()
      if True:
        self.vector_layer = self.pyslip.AddPolygonLayer(
          vector_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<vector_layer>')
        self.vector_text_layer = self.pyslip.AddTextLayer(
          vector_text_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<vector_text_layer>',
          colour='#F62817')
    self.show_filters()

  def get_spotfinder_data(self):
    from scitbx.array_family import flex
    import math
    from dials.algorithms.shoebox import MaskCode
    bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
    fg_code = MaskCode.Valid | MaskCode.Foreground
    strong_code = MaskCode.Valid | MaskCode.Strong

    def map_coords(x, y, p):
      if len(self.pyslip.tiles.raw_image.get_detector()) > 1:
        y, x = self.pyslip.tiles.flex_image.tile_readout_to_picture(
          reflection['panel'], y - 0.5, x - 0.5)
      return self.pyslip.tiles.picture_fast_slow_to_map_relative(
        x, y)

    shoebox_dict = {'width': 2, 'color': '#0000FFA0', 'closed': False}
    ctr_mass_dict = {'width': 2, 'color': '#FF0000', 'closed': False}
    vector_dict = {'width': 4, 'color': '#F62817', 'closed': False}
    i_frame = self.image_chooser.GetClientData(
      self.image_chooser.GetSelection()).index
    imageset = self.image_chooser.GetClientData(
      self.image_chooser.GetSelection()).image_set
    if imageset.get_scan() is not None:
      i_frame += imageset.get_scan().get_array_range()[0]
    shoebox_data = []
    all_pix_data = []
    ctr_mass_data = []
    max_pix_data = []
    predictions_data = []
    miller_indices_data = []
    vector_data = []
    vector_text_data = []
    detector = self.pyslip.tiles.raw_image.get_detector()
    scan = self.pyslip.tiles.raw_image.get_scan()
    to_degrees = 180 / math.pi
    for ref_list in self.reflections:
      if ref_list.has_key('bbox'):
        bbox = ref_list['bbox']
        x0, x1, y0, y1, z0, z1 = bbox.parts()
        bbox_sel = (i_frame >= z0) & (i_frame < z1)
        for reflection in ref_list.select(bbox_sel):
          x0, x1, y0, y1, z0, z1 = reflection['bbox']
          panel = reflection['panel']
          nx = x1 - x0 # size of reflection box in x-direction
          ny = y1 - y0 # size of reflection box in y-direction
          #nz = z1 - z0 # number of frames this spot appears on
          if (self.settings.show_all_pix and reflection.has_key('shoebox')
              and reflection['shoebox'].mask.size() > 0):
            self.show_all_pix_timer.start()
            shoebox = reflection['shoebox']
            iz = i_frame - z0
            print i_frame, z0, iz
            for ix in range(nx):
              for iy in range(ny):
                mask_value = shoebox.mask[iz, iy, ix]
                if ((mask_value == strong_code) or
                    (mask_value == fg_code)):
                  x_, y_ = map_coords(
                    ix + x0 + 0.5, iy + y0 + 0.5, panel)
                  all_pix_data.append((x_, y_))
            self.show_all_pix_timer.stop()

          if self.settings.show_shoebox:
            self.show_shoebox_timer.start()
            x0_, y0_ = map_coords(x0, y0, panel)
            x1_, y1_ = map_coords(x1, y1, panel)
            lines = [(((x0_, y0_), (x0_, y1_)), shoebox_dict),
                     (((x0_, y1_), (x1_, y1_)), shoebox_dict),
                     (((x1_, y1_), (x1_, y0_)), shoebox_dict),
                     (((x1_, y0_), (x0_, y0_)), shoebox_dict)]
            shoebox_data.extend(lines)
            self.show_shoebox_timer.stop()

          if (self.settings.show_max_pix and reflection.has_key('shoebox')
              and reflection['shoebox'].data.size() > 0):
            self.show_max_pix_timer.start()
            shoebox = reflection['shoebox'].data
            offset = flex.max_index(shoebox)
            offset, k = divmod(offset, shoebox.all()[2])
            offset, j = divmod(offset, shoebox.all()[1])
            offset, i = divmod(offset, shoebox.all()[0])
            #assert offset == 0
            max_index = (i, j, k)
            #assert shoebox[max_index] == flex.max(shoebox)
            if z0 + max_index[0] == i_frame:
              x, y = map_coords(x0 + max_index[2] + 0.5,
                                y0 + max_index[1] + 0.5,
                                reflection['panel'])
              max_pix_data.append((x, y))
            self.show_max_pix_timer.stop()

          if (self.settings.show_ctr_mass and
              reflection.has_key('xyzobs.px.value')):
            self.show_ctr_mass_timer.start()
            centroid = reflection['xyzobs.px.value']
            if math.floor(centroid[2]) == i_frame:
              x,y = map_coords(
                centroid[0], centroid[1], reflection['panel'])
              xm1,ym1 = map_coords(
                centroid[0]-1, centroid[1]-1, reflection['panel'])
              xp1,yp1 = map_coords(
                centroid[0]+1, centroid[1]+1, reflection['panel'])
              lines = [(((x, ym1), (x, yp1)), ctr_mass_dict),
                       (((xm1, y), (xp1, y)), ctr_mass_dict)]
              ctr_mass_data.extend(lines)
            self.show_ctr_mass_timer.stop()

      if ((ref_list.has_key('xyzcal.px') or ref_list.has_key('xyzcal.mm')) and
          (self.settings.show_predictions or (
            self.settings.show_miller_indices and ref_list.has_key('miller_index')))):
        if ref_list.has_key('xyzcal.px'):
          frame_numbers = ref_list['xyzcal.px'].parts()[2]
        else:
          phi = ref_list['xyzcal.mm'].parts()[2]
          frame_numbers = scan.get_array_index_from_angle(phi * to_degrees)
        n = 0 # buffer
        frame_predictions_sel = (
          (frame_numbers >= (i_frame-n)) & (frame_numbers < (i_frame+1+n)))
        for reflection in ref_list.select(frame_predictions_sel):
          if (self.settings.show_predictions and
              reflection.has_key('xyzcal.px')):
            x, y = map_coords(reflection['xyzcal.px'][0]+ 0.5,
                              reflection['xyzcal.px'][1] + 0.5,
                              reflection['panel'])
            predictions_data.append((x, y))
          elif (self.settings.show_predictions and
                reflection.has_key('xyzcal.mm')):
            x, y = detector[reflection['panel']].millimeter_to_pixel(
              reflection['xyzcal.mm'][:2])
            x, y = map_coords(x+ 0.5, y + 0.5, reflection['panel'])
            predictions_data.append((x, y))
          if (self.settings.show_miller_indices and
              'miller_index' in reflection and
              reflection['miller_index'] != (0,0,0)):
            miller_indices_data.append((x, y, str(reflection['miller_index']),
                                        {'placement':'ne'}))

    if self.crystals is not None:
      from scitbx import matrix
      from cctbx import crystal
      crystal_model = self.crystals[0]
      cs = crystal.symmetry(unit_cell=crystal_model.get_unit_cell(), space_group=crystal_model.get_space_group())
      cb_op = cs.change_of_basis_op_to_reference_setting()
      crystal_model = crystal_model.change_basis(cb_op)
      A = crystal_model.get_A()
      scan = imageset.get_scan()
      beam = imageset.get_beam()
      phi = scan.get_angle_from_array_index(
        i_frame-imageset.get_array_range()[0], deg=True)
      axis = matrix.col(imageset.get_goniometer().get_rotation_axis())
      beam_centre = detector[0].get_ray_intersection(beam.get_s0())
      beam_x, beam_y = detector[0].millimeter_to_pixel(beam_centre)
      beam_x, beam_y = map_coords(beam_x+ 0.5, beam_y + 0.5, reflection['panel'])
      lines = []
      for i, h in enumerate(((10,0,0), (0,10,0), (0,0,10))):
        r = A * matrix.col(h)
        r_phi = r.rotate_around_origin(axis, phi, deg=True)
        s1 = matrix.col(beam.get_s0()) + r_phi
        x, y = detector[0].get_ray_intersection(s1)
        x, y = detector[0].millimeter_to_pixel((x,y))
        x, y = map_coords(x+ 0.5, y + 0.5, reflection['panel'])
        vector_data.append((((beam_x, beam_y), (x, y)), vector_dict))

        vector_text_data.append((x, y, ('a*', 'b*', 'c*')[i],
                                 {'placement':'ne',
                                  'fontsize': 20,
                                  'color':'#F62817'}))

    from libtbx import group_args
    return group_args(all_pix_data=all_pix_data,
                      shoebox_data=shoebox_data,
                      ctr_mass_data=ctr_mass_data,
                      max_pix_data=max_pix_data,
                      predictions_data=predictions_data,
                      miller_indices_data=miller_indices_data,
                      vector_data=vector_data,
                      vector_text_data=vector_text_data)


class SpotSettingsFrame (SettingsFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SpotSettingsPanel(self, -1)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

class SpotSettingsPanel (SettingsPanel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.settings = self.GetParent().settings
    # CONTROLS 4: additional settings for derived class
    self.settings.show_spotfinder_spots = False
    self.settings.show_dials_spotfinder_spots = True
    self.settings.show_ctr_mass = True
    self.settings.show_max_pix = True
    self.settings.show_all_pix = True
    self.settings.show_shoebox = True
    self.settings.show_predictions = True
    self.settings.show_miller_indices = False
    self.settings.show_mean_filter = False
    self.settings.show_variance_filter = False
    self.settings.show_dispersion = False
    self.settings.show_sigma_b_filter = False
    self.settings.show_sigma_s_filter = False
    self.settings.show_threshold_map = False
    self.settings.nsigma_b = 6
    self.settings.nsigma_s = 3
    self.settings.kernel_size = [3,3]
    self._sizer = wx.BoxSizer(wx.VERTICAL)
    s = self._sizer
    self.SetSizer(self._sizer)
    grid = wx.FlexGridSizer(cols=2, rows=2)
    s.Add(grid)
    txt1 = wx.StaticText(self, -1, "Zoom level:")
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.zoom_ctrl = wx.Choice(self, -1,
      choices=["Auto", "25%", "50%", "100%", "200%", "400%", "800%"])
    self.zoom_ctrl.SetSelection(self.settings.zoom_level)
    grid.Add(self.zoom_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt11 = wx.StaticText(self, -1, "Color scheme:")
    grid.Add(txt11, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.color_ctrl = wx.Choice(self, -1,
      choices=["grayscale","rainbow","heatmap","invert"])
    self.color_ctrl.SetSelection(0)
    grid.Add(self.color_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._sizer.Fit(self)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt2 = wx.StaticText(self, -1, "Brightness")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl = wx.Slider(self, -1, size=(200,-1),
      style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    box.Add(self.brightness_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl.SetMin(1)
    self.brightness_ctrl.SetMax(500)
    self.brightness_ctrl.SetValue(self.settings.brightness)
    self.brightness_ctrl.SetTickFreq(25)

    # Center control
    self.center_ctrl = wx.CheckBox(self, -1, "Mark beam center")
    self.center_ctrl.SetValue(self.settings.show_beam_center)
    s.Add(self.center_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Center of mass control
    self.ctr_mass = wx.CheckBox(self, -1, "Mark centers of mass")
    self.ctr_mass.SetValue(self.settings.show_ctr_mass)
    s.Add(self.ctr_mass, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Max pixel control
    self.max_pix = wx.CheckBox(self, -1, "Spot max pixels")
    self.max_pix.SetValue(self.settings.show_max_pix)
    s.Add(self.max_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot pixels control
    self.all_pix = wx.CheckBox(self, -1, "Spot all pixels")
    self.all_pix.SetValue(self.settings.show_all_pix)
    s.Add(self.all_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot shoebox control
    self.shoebox = wx.CheckBox(self, -1, "Draw reflection shoebox")
    self.shoebox.SetValue(self.settings.show_shoebox)
    s.Add(self.shoebox, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot predictions control
    self.predictions = wx.CheckBox(self, -1, "Show predictions")
    self.predictions.SetValue(self.settings.show_predictions)
    s.Add(self.predictions, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot predictions control
    self.miller_indices = wx.CheckBox(self, -1, "Show hkl")
    self.miller_indices.SetValue(self.settings.show_miller_indices)
    s.Add(self.miller_indices, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Minimum spot area control
    box = wx.BoxSizer(wx.HORIZONTAL)
    from wxtbx.phil_controls.intctrl import IntCtrl
    from wxtbx.phil_controls import EVT_PHIL_CONTROL
    #self.minspotarea_ctrl = IntCtrl(self, -1, pos=(300,180), size=(80,-1),
      #value=self.GetParent().GetParent().horizons_phil.distl.minimum_spot_area,
      #name="Minimum spot area (pxls)")
    #self.minspotarea_ctrl.SetOptional(False)
    #box.Add(self.minspotarea_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    #txtd = wx.StaticText(self, -1,  "Minimum spot area (pxls)",)
    #box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    #s.Add(box)

    # Kabsch thresholding parameters
    grid1 = wx.FlexGridSizer(cols=2, rows=3)
    s.Add(grid1)

    from wxtbx.phil_controls import EVT_PHIL_CONTROL
    from wxtbx.phil_controls.floatctrl import FloatCtrl
    txt1 = wx.StaticText(self, -1, "sigma_background")
    grid1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.nsigma_b_ctrl = FloatCtrl(
      self, value=self.settings.nsigma_b, name="sigma_background")
    self.nsigma_b_ctrl.SetMin(0)
    grid1.Add(self.nsigma_b_ctrl, 0, wx.ALL, 5)

    txt2 = wx.StaticText(self, -1, "sigma_strong")
    grid1.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.nsigma_s_ctrl = FloatCtrl(
      self, value=self.settings.nsigma_s, name="sigma_strong")
    self.nsigma_s_ctrl.SetMin(0)
    grid1.Add(self.nsigma_s_ctrl, 0, wx.ALL, 5)

    from wxtbx.phil_controls.ints import IntsCtrl
    txt3 = wx.StaticText(self, -1, "kernel_size")
    grid1.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.kernel_size_ctrl = IntsCtrl(
      self, value=[3, 3],
      name="kernel_size")
    self.kernel_size_ctrl.SetSize(2)
    self.kernel_size_ctrl.SetMin(1)
    grid1.Add(self.kernel_size_ctrl, 0, wx.ALL, 5)

    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.nsigma_b_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.nsigma_s_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.kernel_size_ctrl)

    from wxtbx.segmentedctrl import SegmentedRadioControl
    self.btn = SegmentedRadioControl(self)
    self.btn.AddSegment("image")
    self.btn.AddSegment("mean")
    self.btn.AddSegment("variance")
    self.btn.AddSegment("dispersion")
    self.btn.AddSegment("sigma_b")
    self.btn.AddSegment("sigma_s ")
    self.btn.AddSegment("threshold")
    self.btn.SetSelection(0)

    self.Bind(wx.EVT_RADIOBUTTON, self.OnUpdateCM, self.btn)
    self.GetSizer().Add(self.btn, 0, wx.ALL, 5)

    btn = wx.Button(self, -1, "Update display", pos=(400, 360))
    self.GetSizer().Add(btn, 0, wx.ALL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnUpdateCM, btn)

    self.collect_values()

    # CONTROLS 3:  Bind events to actions
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.color_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.ctr_mass)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.max_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.all_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.shoebox)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.predictions)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.miller_indices)
    #self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.minspotarea_ctrl)

    txt3 = wx.StaticText(self, -1, "Thumbnail view:")
    s.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.thumb_panel = rstbx.viewer.display.ThumbnailView(
      parent=self,
      size=(256,256),
      style=wx.SUNKEN_BORDER)
    s.Add(self.thumb_panel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  # CONTROLS 2:  Fetch values from widgets
  def collect_values (self) :
    if self.settings.enable_collect_values:
      self.settings.zoom_level = self.zoom_ctrl.GetSelection()
      self.settings.brightness = self.brightness_ctrl.GetValue()
      self.settings.show_beam_center = self.center_ctrl.GetValue()
      self.settings.show_ctr_mass = self.ctr_mass.GetValue()
      self.settings.show_max_pix = self.max_pix.GetValue()
      self.settings.show_all_pix = self.all_pix.GetValue()
      self.settings.show_shoebox = self.shoebox.GetValue()
      self.settings.show_predictions = self.predictions.GetValue()
      self.settings.show_miller_indices = self.miller_indices.GetValue()
      self.settings.color_scheme = self.color_ctrl.GetSelection()
      self.settings.show_mean_filter = self.btn.values[1]
      self.settings.show_variance_filter = self.btn.values[2]
      self.settings.show_dispersion = self.btn.values[3]
      self.settings.show_sigma_b_filter = self.btn.values[4]
      self.settings.show_sigma_s_filter = self.btn.values[5]
      self.settings.show_threshold_map = self.btn.values[6]
      self.settings.nsigma_b = self.nsigma_b_ctrl.GetPhilValue()
      self.settings.nsigma_s = self.nsigma_s_ctrl.GetPhilValue()
      self.settings.kernel_size = self.kernel_size_ctrl.GetPhilValue()

  def OnUpdateCM (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)
