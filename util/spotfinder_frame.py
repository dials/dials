from __future__ import division
import rstbx.viewer.display
import wx
from rstbx.slip_viewer.frame import XrayFrame
from rstbx.viewer.frame import SettingsFrame, SettingsPanel

class SpotFrame(XrayFrame) :
  def __init__ (self, *args, **kwds) :
    self.sweep_filenames = kwds["sweep_filenames"]
    self.reflections = kwds["reflections"]
    del kwds["sweep_filenames"]; del kwds["reflections"] #otherwise wx complains
    super(SpotFrame, self).__init__(*args, **kwds)
    self.viewer.reflections = self.reflections
    self.viewer.frames = self.sweep_filenames
    self.dials_spotfinder_layer = None
    self.shoebox_layer = None
    self.ctr_mass_layer = None

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

  def update_settings(self, layout=True):
    super(SpotFrame, self).update_settings(layout=layout)
    if self.settings.show_dials_spotfinder_spots:
      spotfinder_data = self.get_spotfinder_data()
      shoebox_data = spotfinder_data.shoebox_data
      all_pix_data = spotfinder_data.all_pix_data
      ctr_mass_data = spotfinder_data.ctr_mass_data
      if self.settings.show_all_pix:
        self.dials_spotfinder_layer = self.pyslip.AddPointLayer(
          all_pix_data, color="green", name="<all_pix_layer>",
          radius=2,
          renderer = self.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      elif self.dials_spotfinder_layer is not None:
        self.pyslip.DeleteLayer(self.dials_spotfinder_layer)
        self.dials_spotfinder_layer = None
      if self.settings.show_shoebox:
        self.shoebox_layer = self.pyslip.AddPolygonLayer(
          shoebox_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<shoebox_layer>')
      elif self.shoebox_layer is not None:
        self.pyslip.DeleteLayer(self.shoebox_layer)
        self.shoebox_layer = None
      if self.settings.show_ctr_mass:
        self.ctr_mass_layer = self.pyslip.AddPolygonLayer(
          ctr_mass_data, map_rel=True, visible=True,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<ctr_mass_layer>')
      elif self.ctr_mass_layer is not None:
        self.pyslip.DeleteLayer(self.ctr_mass_layer)
        self.ctr_mass_layer = None

  def get_spotfinder_data(self):

    def map_coords(x, y, p):
      if len(self.pyslip.tiles.raw_image.get_detector()) > 1:
        y, x = self.pyslip.tiles.flex_image.tile_readout_to_picture(
          reflection.panel_number, y - 0.5, x - 0.5)
      return self.pyslip.tiles.picture_fast_slow_to_map_relative(
        x, y)

    shoebox_dict = {'width': 2, 'color': '#0000FFA0', 'closed': False}
    ctr_mass_dict = {'width': 2, 'color': '#FF0000', 'closed': False}
    i_frame = self.image_chooser.GetClientData(
      self.image_chooser.GetSelection()).index
    shoebox_data = []
    all_pix_data = []
    ctr_mass_data = []
    for reflection in self.reflections:
      x0, x1, y0, y1, z0, z1 = reflection.bounding_box
      nx = x1 - x0 # size of reflection box in x-direction
      ny = y1 - y0 # size of reflection box in y-direction
      nz = z1 - z0 # number of frames this spot appears on
      shoebox_mask = reflection.shoebox_mask
      if i_frame >= z0 and i_frame < z1:
        if self.settings.show_all_pix:
          for ix in range(nx):
            for iy in range(ny):
              for iz in range(nz):
                if iz + z0 != i_frame: continue
                if shoebox_mask[iz, iy, ix] > 0:
                  x_, y_ = map_coords(
                    ix + x0 + 0.5, iy + y0 + 0.5, reflection.panel_number)
                  all_pix_data.append((x_, y_))

        if self.settings.show_shoebox:
          x0_, y0_ = map_coords(x0, y0, reflection.panel_number)
          x1_, y1_ = map_coords(x1, y1, reflection.panel_number)
          lines = [(((x0_, y0_), (x0_, y1_)), shoebox_dict),
                   (((x0_, y1_), (x1_, y1_)), shoebox_dict),
                   (((x1_, y1_), (x1_, y0_)), shoebox_dict),
                   (((x1_, y0_), (x0_, y0_)), shoebox_dict)]
          shoebox_data.extend(lines)

        #if self.settings.show_max_pix:
          #pink_brush = wx.Brush(colour=wx.Colour(255,192,203))
          #for spot in self.spotfinder.images[self.frames[0]]["spots_total"]:
            #x,y = self._img.image_coords_as_screen_coords(
              #spot.max_pxl_y(),
              #spot.max_pxl_x())
            #dc.SetPen(wx.Pen('pink'))
            #dc.SetBrush(pink_brush)
            #dc.DrawCircle(x,y,1)

        if self.settings.show_ctr_mass:
          centroid = reflection.centroid_position
          import math
          if math.floor(centroid[2]) == i_frame:
            x,y = map_coords(
              centroid[0], centroid[1], reflection.panel_number)
            xm1,ym1 = map_coords(
              centroid[0]-1, centroid[1]-1, reflection.panel_number)
            xp1,yp1 = map_coords(
              centroid[0]+1, centroid[1]+1, reflection.panel_number)
            lines = [(((x, ym1), (x, yp1)), ctr_mass_dict),
                     (((xm1, y), (xp1, y)), ctr_mass_dict)]
            ctr_mass_data.extend(lines)

    from libtbx import group_args
    return group_args(all_pix_data=all_pix_data,
                      shoebox_data=shoebox_data,
                      ctr_mass_data=ctr_mass_data)


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

    ## Max pixel control
    #self.max_pix = wx.CheckBox(self, -1, "Spot max pixels")
    #self.max_pix.SetValue(self.settings.show_max_pix)
    #s.Add(self.max_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot pixels control
    self.all_pix = wx.CheckBox(self, -1, "Spot all pixels")
    self.all_pix.SetValue(self.settings.show_all_pix)
    s.Add(self.all_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot shoebox control
    self.shoebox = wx.CheckBox(self, -1, "Draw reflection shoebox")
    self.shoebox.SetValue(self.settings.show_shoebox)
    s.Add(self.shoebox, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

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

    self.collect_values()

    # CONTROLS 3:  Bind events to actions
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.color_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.ctr_mass)
    #self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.max_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.all_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.shoebox)
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
      #self.settings.show_max_pix = self.max_pix.GetValue()
      self.settings.show_all_pix = self.all_pix.GetValue()
      self.settings.show_shoebox = self.shoebox.GetValue()
      self.settings.color_scheme = self.color_ctrl.GetSelection()

  def OnUpdateCM (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)
