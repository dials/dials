from __future__ import annotations

import wx


class EllipseSettingsFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
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

        self._pyslip = self.GetParent().GetParent().pyslip
        self._pyslip.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)

        self._points = []
        self._panel = None
        self._point_layer = None

        self.draw_settings()

    def __del__(self):
        try:
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

        title = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(title)
        title_text = wx.StaticText(self, -1, "Click on at least 5 points around a ring")
        title_text.GetFont().SetWeight(wx.BOLD)
        title.Add(title_text, 0, wx.ALL | wx.TOP, border=10)

        grid = wx.FlexGridSizer(cols=4, rows=3, vgap=0, hgap=0)
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
            coords = [
                self._pyslip.tiles.get_flex_pixel_coordinates(*pt)
                for pt in self._points
            ]

            # Filter coords to only those on the same panel
            if len(coords[0]) == 3:
                p = coords[0][2]
                self._panel = int(p)
                coords = [(c[0], c[1]) for c in coords if c[2] == p]

            elif len(coords[0]) == 2:
                self._panel = 0

        if len(coords) >= 5:
            # Dummy values for now
            phi = "0.0"
            l1 = "1.0"
            l2 = "1.0"
            centre = "(0.0, 0.0)"
        else:
            phi = " "
            l1 = " "
            l2 = " "
            centre = " "

        for value in (phi, l1, l2, centre):
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

        # XXX Space for clear and save buttons

        sizer.Layout()
        self.Layout()

    def OnLeftDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPosition()
            self._points.append(self._pyslip.ConvertView2Geo(click_posn))

            # FIXME only draw the point if the last point was on the same panel
            self.DrawPoint()

            self.draw_settings()

        event.Skip()

    def DrawPoint(self):
        if self._point_layer:
            self._pyslip.DeleteLayer(self._point_layer)
            self._point_layer = None

        self._point_layer = self._pyslip.AddPointLayer(
            [(p[0], p[1], {}) for p in self._points],
            name="<predictions_layer>",
            radius=3,
            renderer=self._pyslip.DrawPointLayer,
            color="#00ffff",
            show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
        )
