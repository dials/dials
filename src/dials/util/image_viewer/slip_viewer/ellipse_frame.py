from __future__ import annotations

import numpy as np
import wx
from skimage.measure import EllipseModel

from scitbx import matrix
from wxtbx.phil_controls import EVT_PHIL_CONTROL
from wxtbx.phil_controls.strctrl import StrCtrl


def extract_ellipse_parameters(ellipse: EllipseModel):
    try:
        xc, yc, a, b, theta = ellipse.params
    except AttributeError:
        # Deprecated from skimage 0.26
        (
            xc,
            yc,
        ) = ellipse.center
        a, b = ellipse.axis_lengths
        theta = ellipse.theta

    phi = float(np.degrees(theta))
    a = float(a)
    b = float(b)
    centre_xy = (float(xc), float(yc))

    return phi, a, b, centre_xy


class EllipseSettingsFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.phil_params = args[0].params
        szr = wx.BoxSizer(wx.VERTICAL)
        panel = EllipseSettingsPanel(self)
        self.SetSizer(szr)
        szr.Add(panel, 1, wx.EXPAND)
        szr.Fit(panel)
        self.panel = panel
        self.sizer = szr
        self.Fit()
        self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)

    def Destroy(self):
        for child in self.GetChildren():
            child.Destroy()
        super().Destroy()


class EllipseSettingsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.params = args[0].phil_params

        self._pyslip = self.GetParent().GetParent().pyslip
        self._pyslip.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)

        self._points = []
        self._panel = None
        self._point_layer = None
        self._ellipse_layer = None
        self._show_ellipse = True

        self.draw_settings()

    def __del__(self):
        try:
            self._pyslip.DeleteLayer(self._point_layer)
            self._point_layer = None
            self._pyslip.Unbind(wx.EVT_RIGHT_DOWN, handler=self.OnRightDown)
        except RuntimeError:
            # If the application is closing, the PySlip object has already been deleted
            pass

    def OnClear(self, event):
        self._points = []
        self._pyslip.DeleteLayer(self._point_layer)
        self._point_layer = None
        self._pyslip.DeleteLayer(self._ellipse_layer)
        self._ellipse_layer = None
        self.draw_settings()

    def Destroy(self):
        self.__del__()
        super().Destroy()

    def draw_settings(self):
        for child in self.GetChildren():
            child.Destroy()

        sizer = self.GetSizer()
        if sizer is None:
            sizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(sizer)

        sizer.Clear()

        title = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(title)
        title_text = wx.StaticText(
            self, -1, "Right click on at least 5 points around a ring"
        )
        title_text.GetFont().SetWeight(wx.BOLD)
        title.Add(title_text, 0, wx.ALL | wx.TOP, border=10)

        grid = wx.FlexGridSizer(cols=4, rows=2, vgap=0, hgap=0)
        sizer.Add(grid)

        # Titles
        text = wx.StaticText(self, -1, "phi")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "l1")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "l2")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "centre")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)

        coords = []
        if self._points:
            coords = self._lon_lat_to_slow_fast_panel(self._points)
        self.phi_txt = " "
        self.l1_txt = " "
        self.l2_txt = " "
        self.centre_txt = " "
        enable_save_button = False
        if len(coords) >= 5:
            try:
                ellipse = EllipseModel.from_estimate(np.array(coords))
            except AttributeError:
                # Deprecated from skimage 0.26
                ellipse = EllipseModel()
                success = ellipse.estimate(np.array(coords))
                if not success:
                    ellipse = None
            if not ellipse:
                self.phi_txt = "Fit failed"
            else:
                phi, a, b, centre = extract_ellipse_parameters(ellipse)
                # Use a simplistic model to calculate l1 and l2 scale factors from a and b.
                l1 = 1.0
                l2 = b / a
                self.phi_txt = f"{phi:.2f}"
                self.l1_txt = f"{l1:.6f}"
                self.l2_txt = f"{l2:.6f}"
                self.centre_txt = f"{centre[0]:.2f} {centre[1]:.2f}"
                enable_save_button = True

                # Draw the ellipse
                self._draw_ellipse(ellipse)

        for value in (self.phi_txt, self.l1_txt, self.l2_txt, self.centre_txt):
            grid.Add(
                wx.TextCtrl(
                    self,
                    value=value,
                    size=wx.Size(140, -1),
                    style=wx.TE_READONLY,
                ),
                0,
                wx.ALL,
                5,
            )

        grid = wx.FlexGridSizer(cols=4, rows=1, vgap=0, hgap=0)
        sizer.Add(grid)

        self.clear_button = wx.Button(self, -1, "Clear")
        grid.Add(self.clear_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnClear, self.clear_button)

        self.show_ellipse_ctrl = wx.CheckBox(self, -1, "Show ellipse")
        self.show_ellipse_ctrl.SetValue(self._show_ellipse)
        grid.Add(self.show_ellipse_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 3)
        self.Bind(wx.EVT_CHECKBOX, self.OnShowEllipse, self.show_ellipse_ctrl)

        self.save_params_button = wx.Button(self, -1, "Save parameters")
        self.save_params_button.Enable(enable_save_button)
        grid.Add(self.save_params_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnSaveEllipseParams, self.save_params_button)

        self.save_ellipse_txt_ctrl = StrCtrl(
            self, value=self.params.output.ellipse_params, name="ellipse_phil"
        )
        grid.Add(self.save_ellipse_txt_ctrl, 0, wx.ALL, 5)
        self.Bind(
            EVT_PHIL_CONTROL, self.OnSaveEllipseParams, self.save_ellipse_txt_ctrl
        )

        sizer.Layout()
        self.Layout()

    def _draw_ellipse(self, ellipse: EllipseModel):
        if self._ellipse_layer:
            self._pyslip.DeleteLayer(self._ellipse_layer)
            self._ellipse_layer = None
        phi, a, b, centre = extract_ellipse_parameters(ellipse)
        center = matrix.col(centre)
        e1 = matrix.col((1, 0)).rotate_2d(phi, deg=True)
        e2 = matrix.col((0, 1)).rotate_2d(phi, deg=True)
        ellipse_data = (
            center + a * e1 + b * e2,
            center + a * e1 - b * e2,
            center + a * -e1 - b * e2,
            center + a * -e1 + b * e2,
            center + a * e1 + b * e2,
        )
        ellipse_data = self._slow_fast_to_lon_lat(ellipse_data)

        self._ellipse_layer = self._pyslip.AddEllipseLayer(
            ellipse_data,
            map_rel=True,
            color="#00ffff",
            radius=5,
            visible=True,
            # show_levels=[3,4],
            name="<ellipse_layer>",
        )

    def _lon_lat_to_slow_fast_panel(self, points):
        coords = []
        first_pt = self._pyslip.tiles.get_flex_pixel_coordinates(*points[0])
        if len(first_pt) == 3:
            s, f, self._panel = first_pt
            self._panel = int(self._panel)
        else:
            s, f = coords
            self._panel = 0
        # Correct for half pixel shifts
        coords.append((s + 0.5, f + 0.5))

        for pt in points[1:]:
            coord = self._pyslip.tiles.get_flex_pixel_coordinates(*pt)
            if len(coord) == 3:
                s, f, p = coord
                p = int(p)
            else:
                s, f = coord
                p = 0
            # Skip coordinates not on the same panel as the first point
            if p != self._panel:
                continue
            coords.append((s + 0.5, f + 0.5))

        return coords

    def _slow_fast_to_lon_lat(self, coords):
        points = []
        for coord in coords:
            s, f = coord
            s -= 0.5
            f -= 0.5

            if self._pyslip.tiles.flex_image.supports_rotated_tiles_antialiasing_recommended:
                # Need to undo the effect of picture_to_readout in this case
                s, f = self._pyslip.tiles.flex_image.tile_readout_to_picture(
                    self._panel, s, f
                )

            lon, lat = self._pyslip.tiles.picture_fast_slow_to_map_relative(f, s)

            points.append((lon, lat))
        return points

    def OnRightDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPosition()
            self._points.append(self._pyslip.ConvertView2Geo(click_posn))

            self.DrawPoint()

            self.draw_settings()

        event.Skip()

    def DrawPoint(self):
        if self._point_layer:
            self._pyslip.DeleteLayer(self._point_layer)
            self._point_layer = None

        self._point_layer = self._pyslip.AddPointLayer(
            [(p[0], p[1], {}) for p in self._points],
            name="<points_layer>",
            radius=3,
            renderer=self._pyslip.DrawPointLayer,
            color="#00ffff",
            show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
        )

    def OnShowEllipse(self, event):
        self._show_ellipse = self.show_ellipse_ctrl.GetValue()
        if self._show_ellipse:
            self._pyslip.ShowLayer(self._ellipse_layer)
        else:
            self._pyslip.HideLayer(self._ellipse_layer)

    def OnSaveEllipseParams(self, event):
        self.params.output.ellipse_params = self.save_ellipse_txt_ctrl.GetValue()
        file_name = self.params.output.ellipse_params
        with open(file_name, "w") as f:
            print(f"Saving parameters to {file_name}")
            template = """mode = ellipse
ellipse
{{
  phi = {0}
  l1 = {1}
  l2 = {2}
  centre_xy = {3}
}}
""".format(self.phi_txt, self.l1_txt, self.l2_txt, self.centre_txt)
            f.write(template)
