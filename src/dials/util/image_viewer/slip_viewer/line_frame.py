from __future__ import annotations

import wx


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


class LineSettingsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self._pyslip = self.GetParent().GetParent().pyslip
        self._pyslip.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)

        self._point1 = []
        self._point2 = []
        self._panel = None
        self._line_layer = None

        self.draw_settings()

    def __del__(self):
        self._pyslip.DeleteLayer(self._line_layer)
        self._pyslip.Unbind(wx.EVT_LEFT_DOWN, handler=self.OnLeftDown)

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
        text = wx.StaticText(self, -1, "End")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Mid")
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

        # End
        if self._point2:
            coords = self._pyslip.tiles.get_flex_pixel_coordinates(*self._point2)
            s2, f2 = coords[0:2] + 0.5
            value = f"{f2:.2f},{s2:.2f}"
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

        # Mid
        if self._point1 and self._point2:
            value = f"{(f1 + f2) / 2:.2f},{(s1 + s2) / 2:.2f}"
            # Reset points when the line is finished
            self._point1 = []
            self._point2 = []
            self._panel = None
        else:
            value = ""
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

        sizer.Layout()
        self.Layout()

    def OnLeftDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPosition()

            if not self._point1:
                self._point1 = self._pyslip.ConvertView2Geo(click_posn)
            elif not self._point2:
                self._point2 = self._pyslip.ConvertView2Geo(click_posn)
                self.DrawLine()

            self.draw_settings()

        event.Skip()

    def DrawLine(self):
        if self._line_layer:
            self._pyslip.DeleteLayer(self._line_layer)
            self._line_layer = None

        line_data = [((self._point1, self._point2), {})]

        self._line_layer = self._pyslip.AddPolygonLayer(
            line_data,
            map_rel=True,
            color="#00ffff",
            radius=5,
            visible=True,
            name="<boxsel_pt_layer>",
        )
