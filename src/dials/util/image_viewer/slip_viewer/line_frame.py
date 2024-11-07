from __future__ import annotations

import wx

import wxtbx
from wxtbx import metallicbutton
from wxtbx.phil_controls.strctrl import StrCtrl


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

        self._lines = []
        self._click_points = []
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

        grid = wx.FlexGridSizer(cols=3, vgap=0, hgap=0)
        sizer.Add(grid)
        text = wx.StaticText(self, -1, "Panel:")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        text = wx.StaticText(self, -1, "Line (x1, y1, x2, y2):")
        text.GetFont().SetWeight(wx.BOLD)
        grid.Add(text)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        # Pad lines out to 2 elements for display
        lines = self._lines.copy()
        lines += [(0, None) for i in range(2 - len(lines))]

        for line_id, (panel, points) in enumerate(lines):
            grid.Add(
                StrCtrl(self, value="%i" % (panel), style=wx.TE_READONLY), 0, wx.ALL, 5
            )
            if points:
                value = f"{points[0][0]:.1f} {points[0][1]:.1f} {points[1][0]:.1f} {points[1][1]:.1f}"
            else:
                value = ""
            grid.Add(
                StrCtrl(
                    self,
                    value=value,
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

            self.Bind(
                wx.EVT_BUTTON,
                lambda evt, line_id=line_id: self.OnDeleteLine(evt, line_id=line_id),
                source=btn,
            )

        # self.line_panel_ctrl = IntCtrl(
        #    self, value=0, name="line_panel"
        # )
        # grid.Add(self.line_panel_ctrl, 0, wx.ALL, 5)
        # self.line_ctrl = StrCtrl(self, value="", name="line")
        # grid.Add(self.line_ctrl, 0, wx.ALL, 5)
        # empty cell
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)

        sizer.Layout()
        self.Layout()
        # grid.Layout()
        # sizer.SetSizeHints(self)
        # sizer.Fit(self)

        def OnDeleteLine(self, event, untrusted_id):
            pass

    def OnLeftDown(self, event):
        if not event.ShiftDown():
            click_posn = event.GetPosition()
            xgeo, ygeo = self._pyslip.ConvertView2Geo(click_posn)

            self._click_points.append((xgeo, ygeo))
            if len(self._click_points) > 1:
                self.DrawLine(self._click_points)
                self._lines.append((0, self._click_points))
                self._click_points = []
                self.draw_settings()

        event.Skip()

    def DrawLine(self, vertices):
        if self._line_layer:
            self._pyslip.DeleteLayer(self._line_layer)
            self._line_layer = None

        line_data = []
        d = {}

        for i in range(len(vertices) - 1):
            line_data.append(
                (
                    (vertices[i], vertices[i + 1]),
                    d,
                )
            )

        if line_data:
            self._line_layer = self._pyslip.AddPolygonLayer(
                line_data,
                map_rel=True,
                color="#00ffff",
                radius=5,
                visible=True,
                name="<boxsel_pt_layer>",
            )
