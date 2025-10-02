from __future__ import annotations

import matplotlib
import wx

matplotlib.use("WXAgg")

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

from scitbx import matrix


class LineSettingsFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        szr = wx.BoxSizer(wx.VERTICAL)
        panel = LineSettingsPanel(self)
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


class LineSettingsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self._pyslip = self.GetParent().GetParent().pyslip
        self._pyslip.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)

        self._point1 = []
        self._point2 = []
        self._panel = None
        self._line_layer = None
        self._point_layer = None

        self._lineprof_x = None
        self._lineprof_y = None

        self.draw_settings()

    def __del__(self):
        try:
            self._pyslip.DeleteLayer(self._line_layer)
            self._line_layer = None
            self._pyslip.DeleteLayer(self._point_layer)
            self._point_layer = None
            self._pyslip.Unbind(wx.EVT_LEFT_DOWN, handler=self.OnLeftDown)
        except RuntimeError:
            # If the application is closing, the PySlip object has already been deleted
            pass

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

        grid = wx.FlexGridSizer(cols=4, rows=2, vgap=0, hgap=0)
        sizer.Add(grid)

        # Titles
        text = wx.StaticText(self, -1, "Panel")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Start")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Mid")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "End")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)

        if self._point1:
            coords = self._pyslip.tiles.get_flex_pixel_coordinates(*self._point1)
            if len(coords) == 3:
                s1, f1, p = coords
                self._panel = int(p)
            elif len(coords) == 2:
                s1, f1 = coords
                self._panel = 0
            # Shift to centre of pixel
            s1 += 0.5
            f1 += 0.5

        # Panel
        if self._panel is not None:
            value = f"{self._panel}"
        else:
            value = " "
        grid.Add(
            wx.TextCtrl(self, value=value, size=wx.Size(50, -1), style=wx.TE_READONLY),
            0,
            wx.ALL,
            5,
        )

        # Start
        if self._point1:
            value = f"{f1:.2f},{s1:.2f}"
        else:
            value = " "
        grid.Add(
            wx.TextCtrl(
                self,
                value=value,
                size=wx.Size(130, -1),
                style=wx.TE_READONLY,
            ),
            0,
            wx.ALL,
            5,
        )

        # Mid and end
        if self._point1 and self._point2:
            coords = self._pyslip.tiles.get_flex_pixel_coordinates(*self._point2)
            s2, f2 = coords[0:2] + 0.5
            value_mid = f"{(f1 + f2) / 2:.2f},{(s1 + s2) / 2:.2f}"
            self.calculate_line_profile()
            value_end = f"{f2:.2f},{s2:.2f}"

            # Reset points when the line is finished
            self._point1 = []
            self._point2 = []
            self._panel = None
        else:
            value_mid = " "
            value_end = " "
        grid.Add(
            wx.TextCtrl(
                self,
                value=value_mid,
                size=wx.Size(130, -1),
                style=wx.TE_READONLY,
            ),
            0,
            wx.ALL,
            5,
        )
        grid.Add(
            wx.TextCtrl(
                self,
                value=value_end,
                size=wx.Size(130, -1),
                style=wx.TE_READONLY,
            ),
            0,
            wx.ALL,
            5,
        )

        # Line profile
        figure = Figure(figsize=(1.5, 3))
        figure.subplots_adjust(bottom=0.15)  # Leave room for the x-axis label

        axes = figure.add_subplot(111)
        canvas = FigureCanvas(self, -1, figure)
        sizer.Add(canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        axes.set_title("Line profile")
        axes.set_xlabel("Pixels")
        if self._lineprof_x and self._lineprof_y:
            axes.plot(self._lineprof_x, self._lineprof_y)

        sizer.Layout()
        self.Layout()

    def calculate_line_profile(self):
        pt1 = matrix.col(self._pyslip.tiles.get_flex_pixel_coordinates(*self._point1))
        pt2 = matrix.col(self._pyslip.tiles.get_flex_pixel_coordinates(*self._point2))

        image_data = self._pyslip.tiles.raw_image.get_image_data()
        if not isinstance(image_data, tuple):
            image_data = (image_data,)
        image_data = image_data[self._panel]

        # Sample along the line
        line = pt2 - pt1
        n_samples = int(max(abs(e) for e in line.elems))
        stride = line / n_samples
        step = stride.length()
        x = [step * i for i in range(n_samples + 1)]
        coords = [pt1 + i * stride for i in range(n_samples + 1)]
        vals = []
        for coord in coords:
            isf = (int(coord[0]), int(coord[1]))  # int slow fast
            vals.append(image_data[isf])

        self._lineprof_x = x
        self._lineprof_y = vals

    def OnLeftDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPosition()

            # Update point1 and draw the point
            if not self._point1:
                self._point1 = self._pyslip.ConvertView2Geo(click_posn)
                self.DrawPoint()

            # Or if setting point2, draw the line, if on the same panel
            elif not self._point2:
                pt2 = self._pyslip.ConvertView2Geo(click_posn)
                coords = self._pyslip.tiles.get_flex_pixel_coordinates(*pt2)
                if len(coords) == 3:
                    panel = coords[2]
                else:
                    panel = 0
                if panel == self._panel:
                    self._point2 = pt2
                    self.DrawLine()

            self.draw_settings()

        event.Skip()

    def DrawPoint(self):
        if self._point_layer:
            self._pyslip.DeleteLayer(self._point_layer)
            self._point_layer = None

        self._point_layer = self._pyslip.AddPointLayer(
            [(self._point1[0], self._point1[1], {})],
            name="<predictions_layer>",
            radius=3,
            renderer=self._pyslip.DrawPointLayer,
            color="#00ffff",
            show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
        )

    def DrawLine(self):
        if self._line_layer:
            self._pyslip.DeleteLayer(self._line_layer)
            self._line_layer = None

        self._pyslip.DeleteLayer(self._point_layer)
        self._point_layer = None
        self._pyslip.Update()

        # Calculate bisector position
        pt1 = matrix.col(self._point1)
        pt2 = matrix.col(self._point2)
        mid = (pt1 + pt2) / 2
        v = (pt2 - pt1).normalize()
        orth = v.rotate_2d(90, deg=True)
        pt3 = mid + 10 * orth
        pt4 = mid - 10 * orth

        # Draw line with its bisector
        line_data = [((self._point1, self._point2), {}), ((pt3, pt4), {})]
        self._line_layer = self._pyslip.AddPolygonLayer(
            line_data,
            map_rel=True,
            color="#00ffff",
            radius=5,
            visible=True,
            name="<boxsel_pt_layer>",
        )
