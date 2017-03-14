from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id

import wx

class MaskSettingsFrame(wx.MiniFrame):
  def __init__ (self, *args, **kwds) :
    super(MaskSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.phil_params = args[0].params
    panel = MaskSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class MaskSettingsPanel(wx.Panel):
  def __init__ (self, *args, **kwds) :
    super(MaskSettingsPanel, self).__init__(*args, **kwds)

    self.params = args[0].phil_params

    self._pyslip = self.GetParent().GetParent().pyslip
    self.border_ctrl = None
    self.d_min_ctrl = None
    self.d_max_ctrl = None
    self.draw_settings()
    self.UpdateMask()

  def draw_settings(self):

    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    for child in self.GetChildren():
      if not isinstance(child, FloatSpin):
        # don't destroy FloatSpin controls otherwise bad things happen
        child.Destroy()

    sizer = self.GetSizer()
    if sizer is None:
      sizer = wx.BoxSizer(wx.VERTICAL)
      self.SetSizer(sizer)

    sizer.Clear()
    box = wx.BoxSizer(wx.HORIZONTAL)

    # border control
    if self.border_ctrl is None:
      self.border_ctrl = FloatSpin(
        self, digits=0, name='mask_border', min_val=0)
    box.Add(wx.StaticText(self, label='border'),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.border_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.border_ctrl)
    sizer.Add(box)

    # d_min control
    if self.params.masking.d_min is not None:
      self.d_min = self.params.masking.d_min
    else:
      self.d_min = 0
    box = wx.BoxSizer(wx.HORIZONTAL)
    if self.d_min_ctrl is None:
      self.d_min_ctrl = FloatSpin(
            self, digits=2, name='d_min', value=self.d_min, min_val=0)
    txtd = wx.StaticText(self, label='d_min')
    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.d_min_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.d_min_ctrl)
    sizer.Add(box)

    # d_max control
    if self.params.masking.d_max is not None:
      self.d_max = self.params.masking.d_max
    else:
      self.d_max = 0
    box = wx.BoxSizer(wx.HORIZONTAL)
    if self.d_max_ctrl is None:
      self.d_max_ctrl = FloatSpin(
            self, digits=2, name='d_max', value=self.d_max, min_val=0)
    txtd = wx.StaticText(self, label='d_max')
    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.d_max_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.d_max_ctrl)
    sizer.Add(box)

    untrusted_rectangles = []
    untrusted_polygons = []
    untrusted_circles = []

    for untrusted in self.params.masking.untrusted:
      if untrusted.rectangle is not None:
        untrusted_rectangles.append((untrusted.panel, untrusted.rectangle))
      elif untrusted.polygon is not None:
        untrusted_polygons.append((untrusted.panel, untrusted.polygon))
      elif untrusted.circle is not None:
        untrusted_circles.append((untrusted.panel, untrusted.circle))

    from wxtbx.phil_controls.intctrl import IntCtrl
    from wxtbx.phil_controls.strctrl import StrCtrl
    from wxtbx.phil_controls import EVT_PHIL_CONTROL

    # untrusted rectangles
    grid = wx.FlexGridSizer(cols=2, rows=len(untrusted_rectangles)+2)
    sizer.Add(grid)
    text = wx.StaticText(self, -1, "Panel:")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    text = wx.StaticText(self, -1, "Rectangle (x0, x1, y0, x1):")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    for panel, rectangle in untrusted_rectangles:
      grid.Add(wx.StaticText(self, -1, "%i" %(panel)))
      grid.Add(wx.StaticText(self, -1, "%i %i %i %i" %tuple(rectangle)))
    self.untrusted_rectangle_panel_ctrl = IntCtrl(
      self, value=0, name="untrusted_rectangle_panel")
    grid.Add(self.untrusted_rectangle_panel_ctrl, 0, wx.ALL, 5)
    self.untrusted_rectangle_ctrl = StrCtrl(
      self, value='', name="untrusted_rectangle")
    grid.Add(self.untrusted_rectangle_ctrl, 0, wx.ALL, 5)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_rectangle_ctrl)

    # untrusted polygons
    grid = wx.FlexGridSizer(cols=2, rows=len(untrusted_polygons)+2)
    sizer.Add(grid)
    text = wx.StaticText(self, -1, "Panel:")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    text = wx.StaticText(self, -1, "Polygons (x1, y1, ..., xn, yn):")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    for panel, polygon in untrusted_polygons:
      grid.Add(StrCtrl(self, value="%i" %(panel), style=wx.TE_READONLY), 0, wx.ALL, 5)
      grid.Add(StrCtrl(
        self, value=" ".join(["%i"]*len(polygon)) %tuple(polygon),
        style=wx.TE_READONLY), 0, wx.ALL, 5)
    self.untrusted_polygon_panel_ctrl = IntCtrl(
      self, value=0, name="untrusted_polygon_panel")
    grid.Add(self.untrusted_polygon_panel_ctrl, 0, wx.ALL, 5)
    self.untrusted_polygon_ctrl = StrCtrl(
      self, value='', name="untrusted_polygon")
    grid.Add(self.untrusted_polygon_ctrl, 0, wx.ALL, 5)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_polygon_ctrl)

    # untrusted circles
    grid = wx.FlexGridSizer(cols=2, rows=len(untrusted_circles)+2)
    sizer.Add(grid)
    text = wx.StaticText(self, -1, "Panel:")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    text = wx.StaticText(self, -1, "Circle (x, y, r):")
    text.GetFont().SetWeight(wx.BOLD)
    grid.Add(text)
    for panel, circle in untrusted_circles:
      grid.Add(wx.StaticText(self, -1, "%i" %(panel)))
      grid.Add(wx.StaticText(self, -1, "%i %i %i" %tuple(circle)))
    self.untrusted_circle_panel_ctrl = IntCtrl(
      self, value=0, name="untrusted_circle_panel")
    grid.Add(self.untrusted_circle_panel_ctrl, 0, wx.ALL, 5)
    self.untrusted_circle_ctrl = StrCtrl(
      self, value='', name="untrusted_circle")
    grid.Add(self.untrusted_circle_ctrl, 0, wx.ALL, 5)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_circle_ctrl)

    grid = wx.FlexGridSizer(cols=3, rows=1)
    sizer.Add(grid)

    # show mask control
    self.show_mask_ctrl = wx.CheckBox(self, -1, "Show mask")
    self.show_mask_ctrl.SetValue(self.params.show_mask)
    grid.Add(self.show_mask_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate, self.show_mask_ctrl)

    self.save_mask_button = wx.Button(self, -1, "Save mask")
    grid.Add(self.save_mask_button, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnSaveMask, self.save_mask_button)

    self.save_mask_txt_ctrl = StrCtrl(
      self, value=self.params.output.mask, name="mask_pickle")
    grid.Add(self.save_mask_txt_ctrl, 0, wx.ALL, 5)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.save_mask_txt_ctrl)

    sizer.Layout()
    sizer.Fit(self)

  def __del__(self):
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)

  def OnUpdate(self, event):
    self.params.show_mask = self.show_mask_ctrl.GetValue()
    self.params.output.mask = self.save_mask_txt_ctrl.GetValue()

    if self.d_min_ctrl.GetValue() > 0:
      self.params.masking.d_min = self.d_min_ctrl.GetValue()
    if self.d_max_ctrl.GetValue() > 0:
      self.params.masking.d_max = self.d_max_ctrl.GetValue()
    self.params.masking.border = int(self.border_ctrl.GetValue())

    from dials.util import masking
    from libtbx.utils import flat_list

    untrusted_rectangle = self.untrusted_rectangle_ctrl.GetValue().strip()
    if len(untrusted_rectangle.strip()) > 0:
      rectangle = untrusted_rectangle.strip().replace(',', ' ').split(' ')
      try:
        rectangle = [int(s) for s in rectangle]
        assert len(rectangle) == 4
        panel = int(self.untrusted_rectangle_panel_ctrl.GetValue())
      except Exception, e:
        pass
      else:
        untrusted = masking.phil_scope.extract().untrusted[0]
        untrusted.panel = panel
        untrusted.rectangle = rectangle
        self.params.masking.untrusted.append(untrusted)

    untrusted_polygon = self.untrusted_polygon_ctrl.GetValue().strip()
    if len(untrusted_polygon.strip()) > 0:
      polygon = untrusted_polygon.strip().replace(',', ' ').split(' ')
      try:
        polygon = [int(s) for s in polygon]
        assert len(polygon) % 2 == 0
        assert len(polygon) // 2 > 3
        panel = int(self.untrusted_polygon_panel_ctrl.GetValue())
      except Exception, e:
        pass
      else:
        untrusted = masking.phil_scope.extract().untrusted[0]
        untrusted.panel = panel
        untrusted.polygon = polygon
        self.params.masking.untrusted.append(untrusted)

    untrusted_circle = self.untrusted_circle_ctrl.GetValue().strip()
    if len(untrusted_circle.strip()) > 0:
      circle = untrusted_circle.strip().replace(',', ' ').split(' ')
      try:
        circle = [int(s) for s in circle]
        assert len(circle) == 3
        panel = int(self.untrusted_circle_panel_ctrl.GetValue())
      except Exception, e:
        pass
      else:
        untrusted = masking.phil_scope.extract().untrusted[0]
        untrusted.panel = panel
        untrusted.circle = circle
        self.params.masking.untrusted.append(untrusted)

    self.draw_settings()
    self.UpdateMask()
    image_viewer_frame = self.GetParent().GetParent()
    if image_viewer_frame.settings.show_mask:
      # Force re-drawing of mask
      image_viewer_frame.OnChooseImage(event)

  def OnSaveMask(self, event):
    self.UpdateMask()
    image_viewer_frame = self.GetParent().GetParent()
    mask = image_viewer_frame.mask

    # Save the mask to file
    from libtbx import easy_pickle
    print "Writing mask to %s" % self.params.output.mask
    easy_pickle.dump(self.params.output.mask, mask)

  def UpdateMask(self):

    image_viewer_frame = self.GetParent().GetParent()

    # Generate the mask
    from dials.util.masking import MaskGenerator
    generator = MaskGenerator(self.params.masking)
    imageset = image_viewer_frame.imagesets[0] # XXX
    mask = generator.generate(imageset)

    # Combine with an existing mask, if specified
    if image_viewer_frame.mask is not None:
      for p1, p2 in zip(image_viewer_frame.mask, mask):
        p2 &= p1
      image_viewer_frame.mask = mask
      #self.collect_values()
      image_viewer_frame.update_settings(layout=False)

    image_viewer_frame.mask = mask
