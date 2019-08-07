from __future__ import absolute_import, division, print_function

import wx
from wxtbx.phil_controls.floatctrl import FloatCtrl as _FloatCtrl

# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id

from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin
from wxtbx.phil_controls import EVT_PHIL_CONTROL
from wxtbx.phil_controls.intctrl import IntCtrl
from wxtbx.phil_controls.strctrl import StrCtrl
from wxtbx import metallicbutton
import wxtbx

# Temporary: Make a variable to allow dual API
WX3 = wx.VERSION[0] == 3


class FloatCtrl(_FloatCtrl):

    # override OnFocusLostMethod since calling event.Skip() causes bad things to
    # happen (for reasons I don't understand)
    def OnFocusLost(self, event):
        self.DoSendEvent()
        # event.Skip()


class MaskSettingsFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
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
        self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)


class MaskSettingsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        super(MaskSettingsPanel, self).__init__(*args, **kwds)

        self.params = args[0].phil_params

        self._pyslip = self.GetParent().GetParent().pyslip
        self.border_ctrl = None
        self.d_min_ctrl = None
        self.d_max_ctrl = None
        self.resolution_range_d_min_ctrl = None
        self.resolution_range_d_max_ctrl = None
        self.ice_rings_d_min_ctrl = None
        self.ice_rings_width_ctrl = None
        self._mode_rectangle_layer = None
        self._mode_polygon_layer = None
        self._mode_circle_layer = None
        self._rectangle_x0y0 = None
        self._rectangle_x1y1 = None
        self._mode_rectangle = False
        self._mode_polygon = False
        self._mode_polygon_points = []
        self._mode_circle = False
        self._resolution_range_d_min = 0
        self._resolution_range_d_max = 0

        self._pyslip.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
        self._pyslip.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
        self._pyslip.Bind(wx.EVT_MOTION, self.OnMove)
        self.draw_settings()
        self.UpdateMask()

    def __del__(self):
        if self._mode_rectangle_layer:
            self._pyslip.DeleteLayer(self._mode_rectangle_layer)
        if self._mode_polygon_layer:
            self._pyslip.DeleteLayer(self._mode_polygon_layer)
        if self._mode_circle_layer:
            self._pyslip.DeleteLayer(self._mode_circle_layer)

        self._pyslip.Unbind(wx.EVT_LEFT_DOWN, handler=self.OnLeftDown)
        self._pyslip.Unbind(wx.EVT_LEFT_UP, handler=self.OnLeftUp)
        self._pyslip.Unbind(wx.EVT_MOTION, handler=self.OnMove)

    def draw_settings(self):
        for child in self.GetChildren():
            if not isinstance(child, FloatSpin):
                # don't destroy FloatSpin controls otherwise bad things happen
                child.Destroy()

        sizer = self.GetSizer()
        if sizer is None:
            sizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(sizer)

        sizer.Clear()
        # Put the border, d_min, d_max in an aligned grid
        box = wx.FlexGridSizer(cols=4, hgap=0, vgap=0)

        # border control
        if self.border_ctrl is None:
            self.border_ctrl = FloatSpin(self, digits=0, name="mask_border", min_val=0)
        box.Add(
            wx.StaticText(self, label="Border"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5
        )
        box.Add(
            self.border_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(wx.StaticText(self, label="px"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.border_ctrl)
        # Empty cell after border controls
        box.Add((0, 0))

        # d_min control
        if self.params.masking.d_min is not None:
            self.d_min = self.params.masking.d_min
        else:
            self.d_min = 0
        if self.d_min_ctrl is None:
            self.d_min_ctrl = FloatSpin(
                self,
                digits=2,
                name="d_min",
                value=self.d_min,
                min_val=0,
                increment=0.05,
            )
        txtd = wx.StaticText(self, label="d_min")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(
            self.d_min_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.d_min_ctrl)

        # d_max control
        if self.params.masking.d_max is not None:
            self.d_max = self.params.masking.d_max
        else:
            self.d_max = 0
        if self.d_max_ctrl is None:
            self.d_max_ctrl = FloatSpin(
                self,
                digits=2,
                name="d_max",
                value=self.d_max,
                min_val=0,
                increment=0.05,
            )
        txtd = wx.StaticText(self, label="d_max")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(
            self.d_max_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.d_max_ctrl)
        sizer.Add(box)

        # resolution rings control

        grid = wx.FlexGridSizer(
            cols=2, rows=len(self.params.masking.resolution_range) + 2, vgap=0, hgap=0
        )
        sizer.Add(grid)
        text = wx.StaticText(self, -1, "Resolution range:")
        grid.Add(text)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        for range_id, (d_min, d_max) in enumerate(self.params.masking.resolution_range):
            grid.Add(wx.StaticText(self, -1, "%.2f-%.2f" % (d_min, d_max)))
            btn = metallicbutton.MetallicButton(
                parent=self,
                label="delete",
                bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16),
            )
            grid.Add(btn)
            self.Bind(
                wx.EVT_BUTTON,
                lambda evt, range_id=range_id: self.OnDeleteResolutionRange(
                    evt, range_id=range_id
                ),
                source=btn,
            )

        self.resolution_range_d_min_ctrl = FloatCtrl(
            self, value=self._resolution_range_d_min, name="resolution_range_d_min"
        )
        self.resolution_range_d_min_ctrl.SetMin(0)
        grid.Add(
            self.resolution_range_d_min_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.resolution_range_d_max_ctrl = FloatCtrl(
            self, value=self._resolution_range_d_max, name="resolution_range_d_max"
        )
        self.resolution_range_d_max_ctrl.SetMin(0)
        grid.Add(
            self.resolution_range_d_max_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        # empty cell
        # grid.Add(wx.StaticText(self, -1, ''), 0, wx.EXPAND)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.resolution_range_d_min_ctrl)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.resolution_range_d_max_ctrl)

        # ice rings control
        box = wx.BoxSizer(wx.HORIZONTAL)

        self.ice_rings_ctrl = wx.CheckBox(self, -1, "Ice rings")
        self.ice_rings_ctrl.SetValue(self.params.masking.ice_rings.filter)
        box.Add(self.ice_rings_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_CHECKBOX, self.OnUpdate, self.ice_rings_ctrl)

        # ice rings d_min control
        if self.params.masking.ice_rings.d_min is not None:
            self.ice_rings_d_min = self.params.masking.ice_rings.d_min
        else:
            self.ice_rings_d_min = 0
        if self.ice_rings_d_min_ctrl is None:
            self.ice_rings_d_min_ctrl = FloatSpin(
                self,
                digits=2,
                name="ice_rings_d_min",
                value=self.ice_rings_d_min,
                min_val=0,
                increment=0.05,
            )
        txtd = wx.StaticText(self, label="d_min")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(
            self.ice_rings_d_min_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.ice_rings_d_min_ctrl)

        # ice rings width control
        self.ice_rings_width = self.params.masking.ice_rings.width
        if self.ice_rings_width_ctrl is None:
            self.ice_rings_width_ctrl = FloatSpin(
                self,
                digits=3,
                name="ice_rings_width",
                value=self.ice_rings_width,
                min_val=0.001,
                increment=0.001,
            )
        txtd = wx.StaticText(self, label="width")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(
            self.ice_rings_width_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnUpdate, self.ice_rings_width_ctrl)
        sizer.Add(box)

        untrusted_rectangles = []
        untrusted_polygons = []
        untrusted_circles = []

        # map index in self.params.masking.untrusted to index in above arrays
        self._rectangle_to_untrusted_id = []
        self._polygon_to_untrusted_id = []
        self._circle_to_untrusted_id = []

        for i, untrusted in enumerate(self.params.masking.untrusted):
            if untrusted.rectangle is not None:
                untrusted_rectangles.append((untrusted.panel, untrusted.rectangle))
                self._rectangle_to_untrusted_id.append(i)
            elif untrusted.polygon is not None:
                untrusted_polygons.append((untrusted.panel, untrusted.polygon))
                self._polygon_to_untrusted_id.append(i)
            elif untrusted.circle is not None:
                untrusted_circles.append((untrusted.panel, untrusted.circle))
                self._circle_to_untrusted_id.append(i)

        # untrusted rectangles
        grid = wx.FlexGridSizer(
            cols=3, rows=len(untrusted_rectangles) + 2, vgap=0, hgap=0
        )
        sizer.Add(grid)
        text = wx.StaticText(self, -1, "Panel:")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Rectangle (x0, x1, y0, y1):")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        for rect_id, (panel, rectangle) in enumerate(untrusted_rectangles):
            grid.Add(wx.StaticText(self, -1, "%i" % (panel)))
            grid.Add(wx.StaticText(self, -1, "%i %i %i %i" % tuple(rectangle)))
            btn = metallicbutton.MetallicButton(
                parent=self,
                label="delete",
                bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16),
            )
            grid.Add(btn)
            untrusted_id = self._rectangle_to_untrusted_id[rect_id]
            self.Bind(
                wx.EVT_BUTTON,
                lambda evt, untrusted_id=untrusted_id: self.OnDeleteUntrustedRegion(
                    evt, untrusted_id=untrusted_id
                ),
                source=btn,
            )

        self.untrusted_rectangle_panel_ctrl = IntCtrl(
            self, value=0, name="untrusted_rectangle_panel"
        )
        grid.Add(self.untrusted_rectangle_panel_ctrl, 0, wx.ALL, 5)
        self.untrusted_rectangle_ctrl = StrCtrl(
            self, value="", name="untrusted_rectangle"
        )
        grid.Add(self.untrusted_rectangle_ctrl, 0, wx.ALL, 5)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_rectangle_ctrl)

        # untrusted polygons
        grid = wx.FlexGridSizer(
            cols=3, rows=len(untrusted_polygons) + 2, vgap=0, hgap=0
        )
        sizer.Add(grid)
        text = wx.StaticText(self, -1, "Panel:")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Polygons (x1, y1, ..., xn, yn):")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        for polygon_id, (panel, polygon) in enumerate(untrusted_polygons):
            grid.Add(
                StrCtrl(self, value="%i" % (panel), style=wx.TE_READONLY), 0, wx.ALL, 5
            )
            grid.Add(
                StrCtrl(
                    self,
                    value=" ".join(["%i"] * len(polygon)) % tuple(polygon),
                    style=wx.TE_READONLY,
                ),
                0,
                wx.ALL,
                5,
            )
            btn = metallicbutton.MetallicButton(
                parent=self,
                label="delete",
                bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16),
            )
            grid.Add(btn)
            untrusted_id = self._polygon_to_untrusted_id[polygon_id]
            self.Bind(
                wx.EVT_BUTTON,
                lambda evt, untrusted_id=untrusted_id: self.OnDeleteUntrustedRegion(
                    evt, untrusted_id=untrusted_id
                ),
                source=btn,
            )

        self.untrusted_polygon_panel_ctrl = IntCtrl(
            self, value=0, name="untrusted_polygon_panel"
        )
        grid.Add(self.untrusted_polygon_panel_ctrl, 0, wx.ALL, 5)
        self.untrusted_polygon_ctrl = StrCtrl(self, value="", name="untrusted_polygon")
        grid.Add(self.untrusted_polygon_ctrl, 0, wx.ALL, 5)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_polygon_ctrl)

        # untrusted circles
        grid = wx.FlexGridSizer(cols=3, rows=len(untrusted_circles) + 2, vgap=0, hgap=0)
        sizer.Add(grid)
        text = wx.StaticText(self, -1, "Panel:")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Circle (x, y, r):")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        for circle_id, (panel, circle) in enumerate(untrusted_circles):
            grid.Add(wx.StaticText(self, -1, "%i" % (panel)))
            grid.Add(wx.StaticText(self, -1, "%i %i %i" % tuple(circle)))
            btn = metallicbutton.MetallicButton(
                parent=self,
                label="delete",
                bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16),
            )
            grid.Add(btn)
            untrusted_id = self._circle_to_untrusted_id[circle_id]
            self.Bind(
                wx.EVT_BUTTON,
                lambda evt, untrusted_id=untrusted_id: self.OnDeleteUntrustedRegion(
                    evt, untrusted_id=untrusted_id
                ),
                source=btn,
            )

        self.untrusted_circle_panel_ctrl = IntCtrl(
            self, value=0, name="untrusted_circle_panel"
        )
        grid.Add(self.untrusted_circle_panel_ctrl, 0, wx.ALL, 5)
        self.untrusted_circle_ctrl = StrCtrl(self, value="", name="untrusted_circle")
        grid.Add(self.untrusted_circle_ctrl, 0, wx.ALL, 5)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.untrusted_circle_ctrl)

        # Draw rectangle/circle mode buttons
        grid = wx.FlexGridSizer(cols=4, rows=1, vgap=0, hgap=0)
        sizer.Add(grid)

        grid.Add(wx.StaticText(self, label="Mode:"))
        self.mode_rectangle_button = wx.ToggleButton(self, -1, "Rectangle")
        self.mode_rectangle_button.SetValue(self._mode_rectangle)
        grid.Add(self.mode_rectangle_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnUpdate, self.mode_rectangle_button)

        self.mode_circle_button = wx.ToggleButton(self, -1, "Circle")
        self.mode_circle_button.SetValue(self._mode_circle)
        grid.Add(self.mode_circle_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnUpdate, self.mode_circle_button)

        self.mode_polygon_button = wx.ToggleButton(self, -1, "Polygon")
        self.mode_polygon_button.SetValue(self._mode_polygon)
        grid.Add(self.mode_polygon_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnUpdate, self.mode_polygon_button)

        # show/save mask controls
        grid = wx.FlexGridSizer(cols=3, rows=2, vgap=0, hgap=0)
        sizer.Add(grid)

        self.show_mask_ctrl = wx.CheckBox(self, -1, "Show mask")
        self.show_mask_ctrl.SetValue(self.params.show_mask)
        grid.Add(self.show_mask_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_CHECKBOX, self.OnUpdate, self.show_mask_ctrl)

        self.save_mask_button = wx.Button(self, -1, "Save mask")
        grid.Add(self.save_mask_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnSaveMask, self.save_mask_button)

        self.save_mask_txt_ctrl = StrCtrl(
            self, value=self.params.output.mask, name="mask_pickle"
        )
        grid.Add(self.save_mask_txt_ctrl, 0, wx.ALL, 5)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.save_mask_txt_ctrl)

        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        self.save_params_button = wx.Button(self, -1, "Save")
        grid.Add(self.save_params_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnSaveMaskParams, self.save_params_button)

        self.save_params_txt_ctrl = StrCtrl(
            self, value=self.params.output.mask_params, name="mask_phil"
        )
        grid.Add(self.save_params_txt_ctrl, 0, wx.ALL, 5)
        self.Bind(EVT_PHIL_CONTROL, self.OnUpdate, self.save_params_txt_ctrl)

        sizer.Layout()
        sizer.Fit(self)

    def OnDeleteUntrustedRegion(self, event, untrusted_id):
        if untrusted_id is not None:
            del self.params.masking.untrusted[untrusted_id]
        self.OnUpdate(event)

    def OnDeleteResolutionRange(self, event, range_id):
        if range_id is not None:
            del self.params.masking.resolution_range[range_id]
        self.OnUpdate(event)

    def OnUpdate(self, event):
        image_viewer_frame = self.GetParent().GetParent()

        self.params.show_mask = self.show_mask_ctrl.GetValue()
        # untidy
        image_viewer_frame.settings.show_mask = self.params.show_mask
        image_viewer_frame.params.show_mask = self.params.show_mask
        image_viewer_frame.settings_frame.panel.show_mask.SetValue(
            self.params.show_mask
        )

        self.params.output.mask = self.save_mask_txt_ctrl.GetValue()
        self.params.output.mask_params = self.save_params_txt_ctrl.GetValue()

        if self._mode_polygon and (
            not self.mode_polygon_button.GetValue()
            or self.mode_circle_button.GetValue()
            or self.mode_rectangle_button.GetValue()
        ):
            self.AddUntrustedPolygon(self._mode_polygon_points)
            self._mode_polygon_points = []
            self._pyslip.DeleteLayer(self._mode_polygon_layer)
            self._mode_polygon_layer = None

        if self.mode_rectangle_button.GetValue():
            if not self._mode_rectangle:
                # mode wasn't set but button has been pressed, so set mode
                self._mode_rectangle = True
                # set other modes and buttons to False
                self._mode_circle = False
                self.mode_circle_button.SetValue(False)
                self._mode_polygon = False
                self.mode_polygon_button.SetValue(False)

        if self.mode_circle_button.GetValue():
            if not self._mode_circle:
                # mode wasn't set but button has been pressed, so set mode
                self._mode_circle = True
                # set other modes and buttons to False
                self._mode_rectangle = False
                self.mode_rectangle_button.SetValue(False)
                self._mode_polygon = False
                self.mode_polygon_button.SetValue(False)

        if self.mode_polygon_button.GetValue():
            if not self._mode_polygon:
                # mode wasn't set but button has been pressed, so set mode
                self._mode_polygon = True
                # set other modes and buttons to False
                self._mode_rectangle = False
                self.mode_rectangle_button.SetValue(False)
                self._mode_circle = False
                self.mode_circle_button.SetValue(False)

        if not (
            self.mode_circle_button.GetValue()
            or self.mode_rectangle_button.GetValue()
            or self.mode_polygon_button.GetValue()
        ):
            self._mode_circle = False
            self._mode_rectangle = False
            self._mode_polygon = False

        if self.d_min_ctrl.GetValue() > 0:
            self.params.masking.d_min = self.d_min_ctrl.GetValue()
        else:
            self.params.masking.d_min = None
        if self.d_max_ctrl.GetValue() > 0:
            self.params.masking.d_max = self.d_max_ctrl.GetValue()
        else:
            self.params.masking.d_max = None
        self.params.masking.border = int(self.border_ctrl.GetValue())

        if self.ice_rings_d_min_ctrl.GetValue() > 0:
            self.params.masking.ice_rings.d_min = self.ice_rings_d_min_ctrl.GetValue()
        else:
            self.params.masking.ice_rings.d_min = None
        if self.ice_rings_width_ctrl.GetValue() > 0:
            self.params.masking.ice_rings.width = self.ice_rings_width_ctrl.GetValue()
        self.params.masking.ice_rings.filter = self.ice_rings_ctrl.GetValue()

        self._resolution_range_d_min = float(
            self.resolution_range_d_min_ctrl.GetValue()
        )
        self._resolution_range_d_max = float(
            self.resolution_range_d_max_ctrl.GetValue()
        )

        if self._resolution_range_d_min > 0 and self._resolution_range_d_max > 0:

            self.params.masking.resolution_range.append(
                (self._resolution_range_d_min, self._resolution_range_d_max)
            )
            self._resolution_range_d_min = 0
            self._resolution_range_d_max = 0

        from dials.util import masking
        from libtbx.utils import flat_list

        untrusted_rectangle = self.untrusted_rectangle_ctrl.GetValue().strip()
        if len(untrusted_rectangle.strip()) > 0:
            rectangle = untrusted_rectangle.strip().replace(",", " ").split(" ")
            try:
                rectangle = [int(s) for s in rectangle]
                assert len(rectangle) == 4
                panel = int(self.untrusted_rectangle_panel_ctrl.GetValue())
            except Exception:
                pass
            else:
                untrusted = masking.phil_scope.extract().untrusted[0]
                untrusted.panel = panel
                untrusted.rectangle = rectangle
                self.params.masking.untrusted.append(untrusted)

        untrusted_polygon = self.untrusted_polygon_ctrl.GetValue().strip()
        if len(untrusted_polygon.strip()) > 0:
            polygon = untrusted_polygon.strip().replace(",", " ").split(" ")
            try:
                polygon = [int(s) for s in polygon]
                assert len(polygon) % 2 == 0
                assert len(polygon) // 2 > 3
                panel = int(self.untrusted_polygon_panel_ctrl.GetValue())
            except Exception:
                pass
            else:
                untrusted = masking.phil_scope.extract().untrusted[0]
                untrusted.panel = panel
                untrusted.polygon = polygon
                self.params.masking.untrusted.append(untrusted)

        untrusted_circle = self.untrusted_circle_ctrl.GetValue().strip()
        if len(untrusted_circle.strip()) > 0:
            circle = untrusted_circle.strip().replace(",", " ").split(" ")
            try:
                circle = [int(s) for s in circle]
                assert len(circle) == 3
                panel = int(self.untrusted_circle_panel_ctrl.GetValue())
            except Exception:
                pass
            else:
                untrusted = masking.phil_scope.extract().untrusted[0]
                untrusted.panel = panel
                untrusted.circle = circle
                self.params.masking.untrusted.append(untrusted)

        # Refresh the UI *after* the event is finished
        wx.CallAfter(self.draw_settings)

        self.UpdateMask()

        # Force re-drawing of mask
        image_viewer_frame.OnChooseImage(event)

    def OnSaveMask(self, event):
        self.UpdateMask()
        image_viewer_frame = self.GetParent().GetParent()

        m1 = image_viewer_frame.mask_input
        m2 = image_viewer_frame.mask_image_viewer

        if m1 is not None and m2 is not None:
            mask = []
            for p1, p2 in zip(m1, m2):
                mask.append(p2 & p1)
        elif m1 is not None:
            mask = m1
        elif m2 is not None:
            mask = m2
        else:
            return

        # Save the mask to file
        from libtbx import easy_pickle

        print("Writing mask to %s" % self.params.output.mask)
        easy_pickle.dump(self.params.output.mask, mask)

    def OnSaveMaskParams(self, event):
        from dials.util.masking import phil_scope

        params = phil_scope.extract()
        file_name = self.params.output.mask_params
        with open(file_name, "wb") as f:
            print("Saving parameters to %s" % file_name)
            phil_scope.fetch_diff(phil_scope.format(self.params.masking)).show(f)

    def UpdateMask(self):

        image_viewer_frame = self.GetParent().GetParent()

        # Generate the mask
        from dials.util.masking import MaskGenerator

        generator = MaskGenerator(self.params.masking)
        imageset = image_viewer_frame.imagesets[0]  # XXX
        mask = generator.generate(imageset)

        image_viewer_frame.mask_image_viewer = mask
        image_viewer_frame.update_settings(layout=False)

    def OnLeftDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPositionTuple() if WX3 else event.GetPosition()
            if self._mode_rectangle:
                self._rectangle_x0y0 = click_posn
                self._rectangle_x1y1 = None
                return
            elif self._mode_circle:
                self._circle_xy = click_posn
                self._circle_radius = None
                return
            elif self._mode_polygon:
                self._mode_polygon_points.append(click_posn)
                self.DrawPolygon(self._mode_polygon_points)
        event.Skip()

    def OnLeftUp(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPositionTuple() if WX3 else event.GetPosition()

            if self._mode_rectangle and self._rectangle_x0y0 is not None:
                self._rectangle_x1y1 = click_posn
                x0, y0 = self._rectangle_x0y0
                x1, y1 = self._rectangle_x1y1
                self.AddUntrustedRectangle(x0, y0, x1, y1)
                self._pyslip.DeleteLayer(self._mode_rectangle_layer)
                self._mode_rectangle_layer = None
                self.mode_rectangle_button.SetValue(False)
                self.OnUpdate(event)
                return

            elif self._mode_circle and self._circle_xy is not None:
                xc, yc = self._circle_xy
                xedge, yedge = click_posn
                self.DrawCircle(xc, yc, xedge, yedge)
                try:
                    self.AddUntrustedCircle(xc, yc, xedge, yedge)
                except Exception as e:
                    print(e)
                finally:
                    self._pyslip.DeleteLayer(self._mode_circle_layer)
                    self._mode_circle_layer = None
                    self.mode_circle_button.SetValue(False)
                    self.OnUpdate(event)
                    return

        event.Skip()

    def OnMove(self, event):
        if event.Dragging() and event.LeftIsDown() and not event.ShiftDown():
            if self._mode_rectangle:
                if self._rectangle_x0y0 is not None:
                    x0, y0 = self._rectangle_x0y0
                    x1, y1 = event.GetPositionTuple() if WX3 else event.GetPosition()
                    self.DrawRectangle(x0, y0, x1, y1)
                    return

            elif self._mode_circle:
                if self._circle_xy is not None:
                    xc, yc = self._circle_xy
                    xedge, yedge = (
                        event.GetPositionTuple() if WX3 else event.GetPosition()
                    )
                    self.DrawCircle(xc, yc, xedge, yedge)
                    return
        event.Skip()

    def DrawRectangle(self, x0, y0, x1, y1):
        if self._mode_rectangle_layer:
            self._pyslip.DeleteLayer(self._mode_rectangle_layer)
            self._mode_rectangle_layer = None

        polygon = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
        d = {}
        polygon_data = []
        points = [self._pyslip.ConvertView2Geo(p) for p in polygon]
        for i in range(len(points) - 1):
            polygon_data.append(((points[i], points[i + 1]), d))

        self._mode_rectangle_layer = self._pyslip.AddPolygonLayer(
            polygon_data,
            map_rel=True,
            color="#00ffff",
            radius=5,
            visible=True,
            name="<mode_rectangle_layer>",
        )

    def DrawCircle(self, xc, yc, xedge, yedge):
        if self._mode_circle_layer:
            self._pyslip.DeleteLayer(self._mode_circle_layer)
            self._mode_circle_layer = None

        xc, yc = self._pyslip.ConvertView2Geo((xc, yc))
        xedge, yedge = self._pyslip.ConvertView2Geo((xedge, yedge))

        from scitbx import matrix

        center = matrix.col((xc, yc))
        edge = matrix.col((xedge, yedge))
        r = (center - edge).length()
        if r == 0:
            return

        e1 = matrix.col((1, 0))
        e2 = matrix.col((0, 1))
        circle_data = (
            center + r * (e1 + e2),
            center + r * (e1 - e2),
            center + r * (-e1 - e2),
            center + r * (-e1 + e2),
            center + r * (e1 + e2),
        )

        self._mode_circle_layer = self._pyslip.AddEllipseLayer(
            circle_data,
            map_rel=True,
            color="#00ffff",
            radius=5,
            visible=True,
            # show_levels=[3,4],
            name="<mode_circle_layer>",
        )

    def DrawPolygon(self, vertices):

        if self._mode_polygon_layer:
            self._pyslip.DeleteLayer(self._mode_polygon_layer)
            self._mode_polygon_layer = None

        polygon_data = []
        d = {}

        for i in range(len(vertices) - 1):
            polygon_data.append(
                (
                    (
                        self._pyslip.ConvertView2Geo(vertices[i]),
                        self._pyslip.ConvertView2Geo(vertices[i + 1]),
                    ),
                    d,
                )
            )

        if polygon_data:
            self._mode_polygon_layer = self._pyslip.AddPolygonLayer(
                polygon_data,
                map_rel=True,
                color="#00ffff",
                radius=5,
                visible=True,
                name="<boxsel_pt_layer>",
            )

    def AddUntrustedPolygon(self, vertices):
        if len(vertices) < 4:
            return
        vertices.append(vertices[0])
        vertices = [self._pyslip.ConvertView2Geo(v) for v in vertices]
        vertices = [
            self._pyslip.tiles.map_relative_to_picture_fast_slow(*v) for v in vertices
        ]

        detector = self._pyslip.tiles.raw_image.get_detector()
        if len(detector) > 1:

            point_ = []
            panel_id = None
            for p in vertices:
                p1, p0, p_id = self._pyslip.tiles.flex_image.picture_to_readout(
                    p[1], p[0]
                )
                assert p_id >= 0, "Point must be within a panel"
                if panel_id is not None:
                    assert (
                        panel_id == p_id
                    ), "All points must be contained within a single panel"
                panel_id = p_id
                point_.append((p0, p1))
            vertices = point_

        else:
            panel_id = 0

        from dials.util import masking
        from libtbx.utils import flat_list

        region = masking.phil_scope.extract().untrusted[0]
        points = flat_list(vertices)
        region.polygon = [int(p) for p in points]
        region.panel = panel_id

        self.params.masking.untrusted.append(region)

    def AddUntrustedRectangle(self, x0, y0, x1, y1):
        x0, y0 = self._pyslip.ConvertView2Geo((x0, y0))
        x1, y1 = self._pyslip.ConvertView2Geo((x1, y1))

        if x0 == x1 or y0 == y1:
            return

        points = [(x0, y0), (x1, y1)]
        points = [
            self._pyslip.tiles.map_relative_to_picture_fast_slow(*p) for p in points
        ]

        detector = self._pyslip.tiles.raw_image.get_detector()
        if len(detector) > 1:

            point_ = []
            panel_id = None
            for p in points:
                p1, p0, p_id = self._pyslip.tiles.flex_image.picture_to_readout(
                    p[1], p[0]
                )
                assert p_id >= 0, "Point must be within a panel"
                if panel_id is not None:
                    assert (
                        panel_id == p_id
                    ), "All points must be contained within a single panel"
                panel_id = int(p_id)
                point_.append((p0, p1))
            points = point_

        else:
            panel_id = 0

        (x0, y0), (x1, y1) = points

        if x0 > x1:
            x1, x0 = x0, x1
        if y0 > y1:
            y1, y0 = y0, y1

        panel = detector[panel_id]
        if (
            x1 < 0
            or y1 < 0
            or x0 > panel.get_image_size()[0]
            or y0 > panel.get_image_size()[1]
        ):
            return

        x0 = max(0, x0)
        y0 = max(0, y0)
        x1 = min(panel.get_image_size()[0], x1)
        y1 = min(panel.get_image_size()[1], y1)

        from dials.util import masking
        from libtbx.utils import flat_list

        region = masking.phil_scope.extract().untrusted[0]
        region.rectangle = [int(x0), int(x1), int(y0), int(y1)]
        region.panel = panel_id

        self.params.masking.untrusted.append(region)

    def AddUntrustedCircle(self, xc, yc, xedge, yedge):

        points = [(xc, yc), (xedge, yedge)]

        points = [self._pyslip.ConvertView2Geo(p) for p in points]
        points = [
            self._pyslip.tiles.map_relative_to_picture_fast_slow(*p) for p in points
        ]

        detector = self._pyslip.tiles.raw_image.get_detector()
        if len(detector) > 1:

            points_ = []
            panel_id = None
            for p in points:
                p1, p0, p_id = self._pyslip.tiles.flex_image.picture_to_readout(
                    p[1], p[0]
                )
                assert p_id >= 0, "Point must be within a panel"
                if panel_id is not None:
                    assert (
                        panel_id == p_id
                    ), "All points must be contained within a single panel"
                panel_id = p_id
                points_.append((p0, p1))
            points = points_

        else:
            panel_id = 0

        (xc, yc), (xedge, yedge) = points

        from scitbx import matrix

        center = matrix.col((xc, yc))
        edge = matrix.col((xedge, yedge))
        r = (center - edge).length()
        if r == 0:
            return

        e1 = matrix.col((1, 0))
        e2 = matrix.col((0, 1))
        circle_data = (
            center + r * (e1 + e2),
            center + r * (e1 - e2),
            center + r * (-e1 - e2),
            center + r * (-e1 + e2),
            center + r * (e1 + e2),
        )

        from dials.util import masking
        from libtbx.utils import flat_list

        region = masking.phil_scope.extract().untrusted[0]
        region.circle = [int(xc), int(yc), int(r)]
        region.panel = panel_id

        self.params.masking.untrusted.append(region)
