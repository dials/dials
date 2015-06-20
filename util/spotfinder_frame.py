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
    self._ring_layer = None
    self._resolution_text_layer = None
    self.sel_image_layer = None

    from libtbx.utils import time_log
    self.show_all_pix_timer = time_log("show_all_pix")
    self.show_shoebox_timer = time_log("show_shoebox")
    self.show_max_pix_timer = time_log("show_max_pix")
    self.show_ctr_mass_timer = time_log("show_ctr_mass")
    self.draw_all_pix_timer = time_log("draw_all_pix")
    self.draw_shoebox_timer = time_log("draw_shoebox")
    self.draw_max_pix_timer = time_log("draw_max_pix")
    self.draw_ctr_mass_timer = time_log("draw_ctr_mass_pix")

    self._image_chooser_tmp_key = []
    self._image_chooser_tmp_clientdata = []

    self.init_pyslip_select()

  def init_pyslip_select(self):
    from rstbx.slip_viewer import pyslip
    #self.pyslip.Bind(pyslip.EVT_PYSLIP_SELECT, self.handle_select_event)

    self.TypeMask = 100
    self._xxx_layer = self.pyslip.AddLayer(
      render=self._draw_rings_layer,
      data=[],
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      selectable=True,
      name="<xxx_layer>",
      type=self.TypeMask)
    self.image_layer = self._xxx_layer

    self.add_select_handler(self._xxx_layer, self.boxSelect)
    self.pyslip.SetLayerSelectable(self._xxx_layer, True)

    self.pyslip.layerBSelHandler[self.TypeMask] = self.GetBoxCorners

  def GetBoxCorners(self, layer, p1, p2):
    """Get list of points inside box.

    layer  reference to layer object we are working on
    p1     one corner point of selection box
    p2     opposite corner point of selection box

    We have to figure out which corner is which.

    Return a list of (lon, lat) of points inside box.
    Return None (no selection) or list [((lon, lat), data), ...]
    of points inside the selection box.
    """

    # TODO: speed this up?  Do we need to??
    # get canonical box limits
    (p1x, p1y) = p1
    (p2x, p2y) = p2
    lx = min(p1x, p2x)      # left x coord
    rx = max(p1x, p2x)
    ty = max(p1y, p2y)      # top y coord
    by = min(p1y, p2y)

    return [(lx, by), (lx, ty), (rx, ty), (rx, by)]

  def boxSelect(self, event):
    """Select event from pyslip."""

    from rstbx.slip_viewer import pyslip
    point = event.point
    assert event.evtype == pyslip.EventBoxSelect

    if point:
      assert len(point) == 4

      point = [
        self.pyslip.tiles.map_relative_to_picture_fast_slow(*p)
        for p in point]

      if len(self.pyslip.tiles.raw_image.get_detector()) > 1:
        point = [
          tuple(self.pyslip.tiles.flex_image.picture_to_readout(p[1], p[0]))
          for p in point]
        point = [(p[1],p[0],p[2]) for p in point]
        for p in point: assert p[2] >= 0

      from libtbx.utils import flat_list
      self.settings.untrusted_polygon.append(flat_list(point))

    self.drawUntrustedPolygons()

    return True

  def drawUntrustedPolygons(self):

    # remove any previous selection
    if self.sel_image_layer:
      self.pyslip.DeleteLayer(self.sel_image_layer)
      self.sel_image_layer = None

    untrusted_polygons = self.settings.untrusted_polygon
    data = []
    d = {}
    for polygon in untrusted_polygons:

      if len(self.pyslip.tiles.raw_image.get_detector()) > 1:
        assert len(polygon) % 3 == 0
        polygon = [polygon[i*3:i*3+3] for i in range(len(polygon)//3)]
        polygon = [
          self.pyslip.tiles.flex_image.tile_readout_to_picture(int(p[2]), p[1], p[0])
          for p in polygon]
        polygon = [(p[1], p[0]) for p in polygon]

      else:
        assert len(polygon) % 2 == 0
        polygon = [polygon[i*2:i*2+2] for i in range(len(polygon)//2)]

      points_rel = [self.pyslip.tiles.picture_fast_slow_to_map_relative(*p)
                    for p in polygon]

      points_rel.append(points_rel[0])
      for i in range(len(points_rel)-1):
        data.append(((points_rel[i], points_rel[i+1]), d))

    self.sel_image_layer = \
      self.pyslip.AddPolygonLayer(data, map_rel=True,
                                  color='#00ffff',
                                  radius=5, visible=True,
                                  #show_levels=[3,4],
                                  name='<boxsel_pt_layer>')



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

  def add_file_name_or_data (self, file_name_or_data) :
      """The add_file_name_or_data() function appends @p
      file_name_or_data to the image chooser, unless it is already
      present.  For file-backed images, the base name is displayed in
      the chooser.  If necessary, the number of entries in the chooser
      is pruned.  The function returns the index of the recently added
      entry.  XXX This is probably the place for heuristics to determine
      if the viewer was given a pattern, or a plain list of files.  XXX
      Rename this function, because it only deals with the chooser?
      """

      key = self.get_key(file_name_or_data)
      count = self.image_chooser.GetCount()
      for i in xrange(count) :
        if (key == str(self.image_chooser.GetClientData(i))) :
          return i
      self._image_chooser_tmp_key.append(key)
      self._image_chooser_tmp_clientdata.append(file_name_or_data)
      return len(self._image_chooser_tmp_key) + count

  def load_image (self, file_name_or_data) :
    """The load_image() function displays the image from @p
    file_name_or_data.  The chooser is updated appropriately.
    """

    if len(self._image_chooser_tmp_key):
      starting_count = self.image_chooser.GetCount()
      n = len(self._image_chooser_tmp_key)
      self.image_chooser.AppendItems(self._image_chooser_tmp_key)
      for i in range(n):
        self.image_chooser.SetClientData(
          starting_count+i, self._image_chooser_tmp_clientdata[i])
      self._image_chooser_tmp_key = []
      self._image_chooser_tmp_clientdata = []

    super(SpotFrame, self).load_image(file_name_or_data)

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


  def _draw_rings_layer(self, dc, data, map_rel):
    """Draw a points layer.

    dc       the device context to draw on
    data     an iterable of point tuples:
             (x, y, place, radius, colour, x_off, y_off, pdata)
    map_rel  points relative to map if True, MUST BE TRUE for lightweight
    Assumes all points are the same colour, saving 100's of ms.
    """

    assert map_rel is True
    if len(data)==0:
      return
    (lon, lat, place, radius, colour, x_off, y_off, pdata) = data[0]

    scale = 2**self.pyslip.tiles.zoom_level

    # Draw points on map/view, using transparency if implemented.
    try:
      dc = wx.GCDC(dc)
    except NotImplementedError:
      pass
    dc.SetPen(wx.Pen(colour))
    dc.SetBrush(wx.Brush(colour, wx.TRANSPARENT))
    for (lon, lat, place, radius, colour, x_off, y_off, pdata) in data:
      (x, y) = self.pyslip.ConvertGeo2View((lon, lat))
      dc.DrawCircle(x, y, radius * scale)


  def draw_resolution_rings(self, unit_cell=None, space_group=None):
    from cctbx.array_family import flex
    from cctbx import uctbx
    import math

    image = self.image_chooser.GetClientData(
      self.image_chooser.GetSelection()).image_set
    detector = image.get_detector()
    beam = image.get_beam()

    d_min = detector.get_max_resolution(beam.get_s0())
    d_star_sq_max = uctbx.d_as_d_star_sq(d_min)

    if unit_cell is not None and space_group is not None:
      from cctbx.miller import index_generator
      unit_cell = space_group.average_unit_cell(unit_cell)
      generator = index_generator(unit_cell, space_group.type(), False, d_min)
      indices = generator.to_array()
      spacings = flex.sorted(unit_cell.d(indices))

    else:
      n_rings = 5
      step = d_star_sq_max/n_rings
      spacings = flex.double(
        [uctbx.d_star_sq_as_d((i+1)*step) for i in range(0, n_rings)])

    wavelength = beam.get_wavelength()
    distance = detector[0].get_distance()
    pixel_size = detector[0].get_pixel_size()[0] # FIXME assumes square pixels, and that all panels use same pixel size

    twotheta = uctbx.d_star_sq_as_two_theta(
      uctbx.d_as_d_star_sq(spacings), wavelength)
    L_mm = []
    L_pixels = []
    for tt in twotheta: L_mm.append(distance * math.tan(tt))
    for lmm in L_mm: L_pixels.append(lmm/pixel_size)

    if len(detector) > 1:
      beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(  # FIXME assumes all detector elements use the same
        detector.hierarchy().get_beam_centre(beam.get_s0()))               # millimeter-to-pixel convention
    else:
      beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(
        detector[0].get_beam_centre(beam.get_s0()))

    self._center = [0,0]
    center = self.pyslip.tiles.picture_fast_slow_to_map_relative(
      beam_pixel_fast + self._center[0], beam_pixel_slow + self._center[1])

    # XXX Transparency?
    ring_data = [(center[0], center[1], {"colour": "red", "radius": pxl})
                 for pxl in L_pixels]

    # Remove the old ring layer, and draw a new one.
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self.pyslip.DeleteLayer(self._ring_layer)
      self._ring_layer = None
    self._ring_layer = self.pyslip.AddPointLayer(
      ring_data,
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      renderer=self._draw_rings_layer,
      name="<ring_layer>")

    def rotate_around_point(vector, point, angle, deg=False):
      # http://benn.org/2007/01/06/rotating-coordinates-around-a-centre/
      x, y = vector
      x_centre, y_centre = point
      if deg:
        angle = math.pi * angle / 180
      x_rot = x_centre + math.cos(angle) * (x - x_centre) - math.sin(angle) * (y - y_centre)
      y_rot = y_centre + math.sin(angle) * (x - x_centre) - math.cos(angle) * (y - y_centre)
      return (x_rot, y_rot)

    resolution_text_data = []
    if self.settings.color_scheme > 1 : # heatmap or invert
      textcolour = 'white'
    else:
      textcolour = 'black'
    if unit_cell is None and space_group is None:
      for d, pxl in zip(spacings, L_pixels):
        for angle in (45, 135, 225, 315):
          x, y = rotate_around_point(
            (center[0], center[1]+pxl), center, angle, deg=True)
          resolution_text_data.append(
            (x, y, "%.2f" %d, {
              'placement':'cc', 'colour': 'red',
              'textcolour': textcolour}))

    # Remove the old resolution text layer, and draw a new one.
    if (hasattr(self, "_resolution_text_layer") and self._resolution_text_layer is not None):
      self.pyslip.DeleteLayer(self._resolution_text_layer)
      self._resolution_text_layer = None
    self._resolution_text_layer = self.pyslip.AddTextLayer(
      resolution_text_data,
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      selectable=False,
      name="<resolution_text_layer>",
      colour='red',
      fontsize=15)

  def sum_images(self):
    if self.params.sum_images > 1:
      image = self.pyslip.tiles.raw_image
      detector = image.get_detector()
      if len(detector) == 1:
        raw_data = [image.get_raw_data()]
      else:
        raw_data = [image.get_raw_data(i) for i in range(len(detector))]

      i_frame = self.image_chooser.GetClientData(
        self.image_chooser.GetSelection()).index
      imageset = self.image_chooser.GetClientData(i_frame).image_set

      for i in range(1, self.params.sum_images):
        if (i_frame + i) >= len(imageset): break
        raw_data_i = [imageset[i_frame + i]]
        for j, rd in enumerate(raw_data):
          rd += raw_data_i[j]

      self.pyslip.tiles.set_image_data(raw_data)

      self.pyslip.ZoomToLevel(self.pyslip.tiles.zoom_level)
      self.update_statusbar() # XXX Not always working?
      self.Layout()

  def show_filters(self):
    from dials.algorithms.image.threshold import KabschDebug
    from dials.array_family import flex

    image = self.pyslip.tiles.raw_image
    detector = image.get_detector()
    if len(detector) == 1:
      raw_data = [image.get_raw_data()]
    else:
      raw_data = [image.get_raw_data(i) for i in range(len(detector))]

    if (self.settings.show_mean_filter or
        self.settings.show_variance_filter or
        self.settings.show_dispersion or
        self.settings.show_sigma_b_filter or
        self.settings.show_sigma_s_filter or
        self.settings.show_global_threshold_filter or
        self.settings.show_threshold_map):

      trange = [p.get_trusted_range() for p in detector]
      mask = []
      for tr, im in zip(trange, raw_data):
        mask.append(im > int(tr[0]))

      gain_value = self.settings.gain
      assert gain_value > 0
      gain_map = [flex.double(raw_data[i].accessor(), gain_value)
                  for i in range(len(detector))]

      nsigma_b = self.settings.nsigma_b
      nsigma_s = self.settings.nsigma_s
      global_threshold = self.settings.global_threshold
      min_local = self.settings.min_local
      size = self.settings.kernel_size
      kabsch_debug_list = []
      for i_panel in range(len(detector)):
        kabsch_debug_list.append(
          KabschDebug(
            raw_data[i_panel].as_double(), mask[i_panel], gain_map[i_panel],
            size, nsigma_b, nsigma_s, global_threshold, min_local))

      if self.settings.show_mean_filter:
        mean = [kabsch.mean() for kabsch in kabsch_debug_list]
        self.pyslip.tiles.set_image_data(mean)
      elif self.settings.show_variance_filter:
        variance = [kabsch.variance() for kabsch in kabsch_debug_list]
        self.pyslip.tiles.set_image_data(variance)
      elif self.settings.show_dispersion:
        cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
        self.pyslip.tiles.set_image_data(cv)
      elif self.settings.show_sigma_b_filter:
        cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
        cv_mask = [kabsch.cv_mask() for kabsch in kabsch_debug_list]
        cv_mask = [mask.as_1d().as_double() for mask in cv_mask]
        for i, mask in enumerate(cv_mask):
          mask.reshape(cv[i].accessor())
        self.pyslip.tiles.set_image_data(cv_mask)
      elif self.settings.show_sigma_s_filter:
        cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
        value_mask = [kabsch.value_mask() for kabsch in kabsch_debug_list]
        value_mask = [mask.as_1d().as_double() for mask in value_mask]
        for i, mask in enumerate(value_mask):
          mask.reshape(cv[i].accessor())
        self.pyslip.tiles.set_image_data(value_mask)
      elif self.settings.show_global_threshold_filter:
        cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
        global_mask = [kabsch.global_mask() for kabsch in kabsch_debug_list]
        global_mask = [mask.as_1d().as_double() for mask in global_mask]
        for i, mask in enumerate(global_mask):
          mask.reshape(cv[i].accessor())
        self.pyslip.tiles.set_image_data(global_mask)
      elif self.settings.show_threshold_map:
        cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
        final_mask = [kabsch.final_mask() for kabsch in kabsch_debug_list]
        final_mask = [mask.as_1d().as_double() for mask in final_mask]
        for i, mask in enumerate(final_mask):
          mask.reshape(cv[i].accessor())
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
      if self._ring_layer is not None:
        self.pyslip.DeleteLayer(self._ring_layer)
        self._ring_layer = None
      if self._resolution_text_layer is not None:
        self.pyslip.DeleteLayer(self._resolution_text_layer)
        self._resolution_text_layer = None

      if self.settings.show_miller_indices and len(miller_indices_data):
        self.miller_indices_layer = self.pyslip.AddTextLayer(
          miller_indices_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<miller_indices_layer>')
      if self.settings.show_predictions and len(predictions_data):
        self.predictions_layer = self.pyslip.AddPointLayer(
          predictions_data, name="<predictions_layer>",
          radius=3,
          renderer = self.pyslip.DrawPointLayer,
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
    self.sum_images()
    if self.params.sum_images == 1:
      self.show_filters()
    if self.settings.show_resolution_rings:
      self.draw_resolution_rings()
    elif self.settings.show_ice_rings:
      from cctbx import sgtbx, uctbx
      unit_cell = uctbx.unit_cell((4.498,4.498,7.338,90,90,120))
      space_group = sgtbx.space_group_info(number=194).group()
      self.draw_resolution_rings(unit_cell=unit_cell, space_group=space_group)
    self.drawUntrustedPolygons()

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
    #prediction_colours = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                          #"#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                          #"#cab2d6"] * 10
    # alternative colour scheme
    prediction_colours = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                          "#ff7f00", "#ffff33", "#a65628", "#f781bf",
                          "#999999"] * 10
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
            #print i_frame, z0, iz
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
        for i_expt in range(flex.max(ref_list['id'])+1):
          expt_sel = ref_list['id'] == i_expt
          frame_predictions_sel = (
            (frame_numbers >= (i_frame-n)) & (frame_numbers < (i_frame+1+n)))
          for reflection in ref_list.select(frame_predictions_sel & expt_sel):
            if (self.settings.show_predictions and
                reflection.has_key('xyzcal.px')):
              x, y = map_coords(reflection['xyzcal.px'][0]+ 0.5,
                                reflection['xyzcal.px'][1] + 0.5,
                                reflection['panel'])
              predictions_data.append(
                (x, y, {'colour':prediction_colours[i_expt]}))
            elif (self.settings.show_predictions and
                  reflection.has_key('xyzcal.mm')):
              x, y = detector[reflection['panel']].millimeter_to_pixel(
                reflection['xyzcal.mm'][:2])
              x, y = map_coords(x+ 0.5, y + 0.5, reflection['panel'])
              predictions_data.append(
                (x, y, {'colour':prediction_colours[i_expt]}))
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
      beam_x, beam_y = map_coords(beam_x+ 0.5, beam_y + 0.5, 0)
      lines = []
      for i, h in enumerate(((10,0,0), (0,10,0), (0,0,10))):
        r = A * matrix.col(h)
        r_phi = r.rotate_around_origin(axis, phi, deg=True)
        s1 = matrix.col(beam.get_s0()) + r_phi
        x, y = detector[0].get_ray_intersection(s1)
        x, y = detector[0].millimeter_to_pixel((x,y))
        x, y = map_coords(x+ 0.5, y + 0.5, 0)
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

  def get_detector(self):
    return self.imagesets[0].get_detector()

  def get_beam(self):
    return self.imagesets[0].get_beam()


class SpotSettingsFrame (SettingsFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    self.params = self.GetParent().params
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
    self.params = self.GetParent().params
    # CONTROLS 4: additional settings for derived class
    self.settings.brightness = self.params.brightness
    self.settings.color_scheme = self.params.color_scheme
    self.settings.show_spotfinder_spots = False
    self.settings.show_dials_spotfinder_spots = True
    self.settings.show_resolution_rings = self.params.show_resolution_rings
    self.settings.untrusted_polygon = self.params.untrusted_polygon
    self.settings.show_ice_rings = self.params.show_ice_rings
    self.settings.show_ctr_mass = self.params.show_ctr_mass
    self.settings.show_max_pix = self.params.show_max_pix
    self.settings.show_all_pix = self.params.show_all_pix
    self.settings.show_shoebox = self.params.show_shoebox
    self.settings.show_predictions = self.params.show_predictions
    self.settings.show_miller_indices = self.params.show_miller_indices
    self.settings.show_mean_filter = self.params.display == "mean"
    self.settings.show_variance_filter = self.params.display == "variance"
    self.settings.show_dispersion = self.params.display == "dispersion"
    self.settings.show_sigma_b_filter = self.params.display == "sigma_b"
    self.settings.show_sigma_s_filter = self.params.display == "sigma_s"
    self.settings.show_threshold_map = self.params.display == "threshold"
    self.settings.show_global_threshold_filter = self.params.display == "global_threshold"
    self.settings.nsigma_b = self.params.nsigma_b
    self.settings.nsigma_s = self.params.nsigma_s
    self.settings.global_threshold = self.params.global_threshold
    self.settings.kernel_size = self.params.kernel_size
    self.settings.min_local = self.params.min_local
    self.settings.gain = self.params.gain
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
    color_schemes = ["grayscale","rainbow","heatmap","invert"]
    self.color_ctrl = wx.Choice(self, -1,
      choices=color_schemes)
    self.color_ctrl.SetSelection(color_schemes.index(self.params.color_scheme))
    grid.Add(self.color_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._sizer.Fit(self)

    from wxtbx.phil_controls.intctrl import IntCtrl
    from wxtbx.phil_controls import EVT_PHIL_CONTROL

    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    grid = wx.FlexGridSizer(cols=1, rows=2)
    box.Add(grid)
    txt2 = wx.StaticText(self, -1, "Brightness:")
    grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_txt_ctrl = IntCtrl(
      self, value=self.settings.brightness, name="brightness")
    grid.Add(self.brightness_txt_ctrl, 0, wx.ALL, 5)
    self.brightness_txt_ctrl.SetMin(1)
    self.brightness_txt_ctrl.SetMax(500)
    self.brightness_ctrl = wx.Slider(self, -1, size=(150,-1),
      style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    self.brightness_ctrl.SetMin(1)
    self.brightness_ctrl.SetMax(500)
    self.brightness_ctrl.SetValue(self.settings.brightness)
    self.brightness_ctrl.SetTickFreq(25)
    box.Add(self.brightness_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    grid = wx.FlexGridSizer(cols=2, rows=8)
    s.Add(grid)

    # Resolution rings control
    self.resolution_rings_ctrl = wx.CheckBox(self, -1, "Show resolution rings")
    self.resolution_rings_ctrl.SetValue(self.settings.show_resolution_rings)
    grid.Add(self.resolution_rings_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Ice rings control
    self.ice_rings_ctrl = wx.CheckBox(self, -1, "Show ice rings")
    self.ice_rings_ctrl.SetValue(self.settings.show_ice_rings)
    grid.Add(self.ice_rings_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Center control
    self.center_ctrl = wx.CheckBox(self, -1, "Mark beam center")
    self.center_ctrl.SetValue(self.settings.show_beam_center)
    grid.Add(self.center_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Center of mass control
    self.ctr_mass = wx.CheckBox(self, -1, "Mark centers of mass")
    self.ctr_mass.SetValue(self.settings.show_ctr_mass)
    grid.Add(self.ctr_mass, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Max pixel control
    self.max_pix = wx.CheckBox(self, -1, "Spot max pixels")
    self.max_pix.SetValue(self.settings.show_max_pix)
    grid.Add(self.max_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot pixels control
    self.all_pix = wx.CheckBox(self, -1, "Spot all pixels")
    self.all_pix.SetValue(self.settings.show_all_pix)
    grid.Add(self.all_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot shoebox control
    self.shoebox = wx.CheckBox(self, -1, "Draw reflection shoebox")
    self.shoebox.SetValue(self.settings.show_shoebox)
    grid.Add(self.shoebox, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot predictions control
    self.predictions = wx.CheckBox(self, -1, "Show predictions")
    self.predictions.SetValue(self.settings.show_predictions)
    grid.Add(self.predictions, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot predictions control
    self.miller_indices = wx.CheckBox(self, -1, "Show hkl")
    self.miller_indices.SetValue(self.settings.show_miller_indices)
    grid.Add(self.miller_indices, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    self.clear_all_button = wx.Button(self, -1, "Clear all")
    grid.Add(self.clear_all_button, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnClearAll, self.clear_all_button)

    self.save_mask_button = wx.Button(self, -1, "Save mask")
    grid.Add(self.save_mask_button, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnSaveMask, self.save_mask_button)

    # Minimum spot area control
    box = wx.BoxSizer(wx.HORIZONTAL)
    #self.minspotarea_ctrl = IntCtrl(self, -1, pos=(300,180), size=(80,-1),
      #value=self.GetParent().GetParent().horizons_phil.distl.minimum_spot_area,
      #name="Minimum spot area (pxls)")
    #self.minspotarea_ctrl.SetOptional(False)
    #box.Add(self.minspotarea_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    #txtd = wx.StaticText(self, -1,  "Minimum spot area (pxls)",)
    #box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    #s.Add(box)

    # Kabsch thresholding parameters
    grid1 = wx.FlexGridSizer(cols=2, rows=6)
    s.Add(grid1)

    from wxtbx.phil_controls.floatctrl import FloatCtrl
    txt1 = wx.StaticText(self, -1, "Sigma background")
    grid1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.nsigma_b_ctrl = FloatCtrl(
      self, value=self.settings.nsigma_b, name="sigma_background")
    self.nsigma_b_ctrl.SetMin(0)
    grid1.Add(self.nsigma_b_ctrl, 0, wx.ALL, 5)

    txt2 = wx.StaticText(self, -1, "Sigma strong")
    grid1.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.nsigma_s_ctrl = FloatCtrl(
      self, value=self.settings.nsigma_s, name="sigma_strong")
    self.nsigma_s_ctrl.SetMin(0)
    grid1.Add(self.nsigma_s_ctrl, 0, wx.ALL, 5)

    txt1 = wx.StaticText(self, -1, "Global Threshold")
    grid1.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.global_threshold_ctrl = FloatCtrl(
      self, value=self.settings.global_threshold, name="global_threshold")
    self.global_threshold_ctrl.SetMin(0)
    grid1.Add(self.global_threshold_ctrl, 0, wx.ALL, 5)

    from wxtbx.phil_controls.intctrl import IntCtrl
    txt4 = wx.StaticText(self, -1, "Min. local")
    grid1.Add(txt4, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.min_local_ctrl = IntCtrl(
      self, value=self.settings.min_local, name="min_local")
    self.min_local_ctrl.SetMin(0)
    grid1.Add(self.min_local_ctrl, 0, wx.ALL, 5)

    txt4 = wx.StaticText(self, -1, "Gain")
    grid1.Add(txt4, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.gain_ctrl = FloatCtrl(
      self, value=self.settings.gain, name="gain")
    self.gain_ctrl.SetMin(0)
    grid1.Add(self.gain_ctrl, 0, wx.ALL, 5)

    from wxtbx.phil_controls.ints import IntsCtrl
    txt3 = wx.StaticText(self, -1, "Kernel size")
    grid1.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.kernel_size_ctrl = IntsCtrl(
      self, value=[3, 3],
      name="kernel_size")
    self.kernel_size_ctrl.SetSize(2)
    self.kernel_size_ctrl.SetMin(1)
    grid1.Add(self.kernel_size_ctrl, 0, wx.ALL, 5)

    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.brightness_txt_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.nsigma_b_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.nsigma_s_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.global_threshold_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.kernel_size_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.min_local_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.gain_ctrl)

    from wxtbx.segmentedctrl import SegmentedRadioControl, SEGBTN_VERTICAL
    self.btn = SegmentedRadioControl(self, style=SEGBTN_VERTICAL)
    self.btn.AddSegment("image")
    self.btn.AddSegment("mean")
    self.btn.AddSegment("variance")
    self.btn.AddSegment("dispersion")
    self.btn.AddSegment("global")
    self.btn.AddSegment("sigma_b")
    self.btn.AddSegment("sigma_s ")
    self.btn.AddSegment("threshold")
    self.btn.SetSelection(
      ["image", "mean", "variance", "dispersion", "global_threshold",
       "sigma_b", "sigma_s", "threshold"].index(self.params.display))

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
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.resolution_rings_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.ice_rings_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.ctr_mass)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.max_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.all_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.shoebox)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.predictions)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.miller_indices)
    #self.Bind(EVT_PHIL_CONTROL, self.OnUpdateCM, self.minspotarea_ctrl)

    have_thumbnail = False
    if have_thumbnail:
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
      self.settings.show_resolution_rings = self.resolution_rings_ctrl.GetValue()
      self.settings.show_ice_rings = self.ice_rings_ctrl.GetValue()
      self.settings.zoom_level = self.zoom_ctrl.GetSelection()
      # get brightness from slider or text box, whichever is different from
      # the current value, then update the other input field so that both
      # display the current value
      if self.brightness_ctrl.GetValue() != self.settings.brightness:
        self.settings.brightness = self.brightness_ctrl.GetValue()
        self.brightness_txt_ctrl.SetValue(self.settings.brightness)
      else:
        self.settings.brightness = int(self.brightness_txt_ctrl.GetValue())
        self.brightness_ctrl.SetValue(self.settings.brightness)
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
      self.settings.show_global_threshold_filter = self.btn.values[4]
      self.settings.show_sigma_b_filter = self.btn.values[5]
      self.settings.show_sigma_s_filter = self.btn.values[6]
      self.settings.show_threshold_map = self.btn.values[7]
      self.settings.nsigma_b = self.nsigma_b_ctrl.GetPhilValue()
      self.settings.nsigma_s = self.nsigma_s_ctrl.GetPhilValue()
      self.settings.global_threshold = self.global_threshold_ctrl.GetPhilValue()
      self.settings.kernel_size = self.kernel_size_ctrl.GetPhilValue()
      self.settings.min_local = self.min_local_ctrl.GetPhilValue()
      self.settings.gain = self.gain_ctrl.GetPhilValue()

  def OnUpdateCM (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)

  def OnClearAll(self, event):
    for btn in (self.center_ctrl, self.ctr_mass, self.max_pix, self.all_pix,
                self.shoebox, self.predictions, self.miller_indices,
                self.ice_rings_ctrl, self.resolution_rings_ctrl):
      btn.SetValue(False)
    self.OnUpdateCM(event)

  def OnSaveMask(self, event):
    print "Saving mask"
    from scitbx.array_family import flex

    imagesets = self.GetParent().GetParent().imagesets # XXX
    detector = imagesets[0].get_detector()

    from dials.algorithms.peak_finding.spotfinder_factory import polygon
    polygons = self.settings.untrusted_polygon

    if len(detector) > 1:
      polygons = [[vertices[i*3:i*3+3] for i in range(len(vertices)//3)] for vertices in polygons]
      panel_ids = [poly[0][2] for poly in polygons]
      polygons = [polygon([p[:2] for p in poly]) for poly in polygons]

    else:
      panel_ids = None
      polygons = [
        polygon([vertices[i*2:i*2+2] for i in range(len(vertices)//2)]) for vertices in polygons]

    # Create the mask for each image
    masks = []
    # Get the first image
    image = imagesets[0][0]
    if not isinstance(image, tuple):
      image = (image,)

    for i_panel, (im, panel) in enumerate(zip(image, detector)):
      mask = flex.bool(flex.grid(im.all()), True)

      import math
      for i, poly in enumerate(polygons):
        if panel_ids is not None and panel_ids[i] != i_panel:
          continue
        min_x = int(math.floor(min(v[0] for v in poly.vertices)))
        max_x = int(math.ceil(max(v[0] for v in poly.vertices)))
        min_y = int(math.floor(min(v[1] for v in poly.vertices)))
        max_y = int(math.ceil(max(v[1] for v in poly.vertices)))

        for i in range(min_x, max_x+1):
          for j in range(min_y, max_y+1):
            if poly.is_inside(i,j):
              mask[j,i] = False

      # Add to the list
      masks.append(mask)
      #print mask.count(True), mask.count(False)

    from libtbx import easy_pickle
    easy_pickle.dump('mask.pickle', tuple(masks))

    print "Saved mask.pickle"
