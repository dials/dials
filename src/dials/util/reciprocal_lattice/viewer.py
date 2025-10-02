from __future__ import annotations

from math import pi

import wx
from annlib_ext import AnnAdaptorSelfInclude
from wx.lib.agw import floatspin

import gltbx
import gltbx.gl as gl
import libtbx.phil
import wxtbx.utils
from libtbx import Auto
from scitbx.array_family import flex
from scitbx.math import minimum_covering_sphere
from wxtbx.segmentedctrl import (
    SEGBTN_HORIZONTAL,
    SegmentedRadioControl,
    SegmentedToggleControl,
)

from dials.util import wx_viewer
from dials.util.reciprocal_lattice import Render3d

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.reciprocal_lattice.phil_scope

show_rotation_axis = False
  .type = bool
show_beam_vector = False
  .type = bool
show_reciprocal_cell = False
  .type = bool
label_nearest_point = False
  .type = bool
marker_size = Auto
  .type = int(value_min=1)
autospin = False
  .type = bool
model_view_matrix = None
  .type = floats(size=16)
max_experiments = 15
  .type = int(value_min=0)
""",
    process_includes=True,
)


# WX3 - WX4 compatibility
def _rewrite_event(unbound):
    """Decorator to intercept the event and add missing instance methods"""

    def _wrapp(self, event):
        event.GetPositionTuple = event.GetPosition
        return unbound(self, event)

    return _wrapp


# HACK: Monkeypatch wxtbx so that we don't use old interfaces
wxtbx.segmentedctrl.SegmentedControl.HitTest = _rewrite_event(
    wxtbx.segmentedctrl.SegmentedControl.HitTest
)
wxtbx.segmentedctrl.SegmentedControl.OnMotion = _rewrite_event(
    wxtbx.segmentedctrl.SegmentedControl.OnMotion
)


class ReciprocalLatticeViewer(wx.Frame, Render3d):
    def __init__(self, parent, id, title, size, settings=None, *args, **kwds):
        wx.Frame.__init__(self, parent, id, title, size=size, *args, **kwds)
        Render3d.__init__(self, settings=settings)
        self.parent = self.GetParent()
        self.statusbar = self.CreateStatusBar()
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.create_settings_panel()
        self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
        self.create_viewer_panel()
        self.sizer.Add(self.viewer, 1, wx.EXPAND | wx.ALL)
        self.SetSizerAndFit(self.sizer)
        self.SetMinSize(self.settings_panel.GetSize())
        self.Bind(wx.EVT_CLOSE, self.OnClose, self)
        self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
        self.Bind(wx.EVT_ACTIVATE, self.OnActive)
        self.viewer.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
        self.viewer.SetFocus()

    def load_models(self, experiments, reflections):
        Render3d.load_models(self, experiments, reflections)
        if self.settings.beam_centre is not None:
            self.settings_panel.beam_panel_ctrl.SetValue(
                self.settings.beam_centre_panel
            )
            self.settings_panel.beam_fast_ctrl.SetValue(self.settings.beam_centre[0])
            self.settings_panel.beam_slow_ctrl.SetValue(self.settings.beam_centre[1])
        else:
            self.settings_panel.beam_fast_ctrl.Disable()
            self.settings_panel.beam_slow_ctrl.Disable()
        if self.settings.marker_size is Auto:
            max_radius = max(self.reflections["rlp"].norms())
            volume = 4 / 3 * pi * max_radius**3
            density = len(self.reflections) / volume
            # Set marker size to between 5 and 50 depending on density, where
            # 1000 < density < 20000 ==> 50 < marker_size < 5
            marker_size = (-45 / 19000) * density + (5 + 900 / 19)
            marker_size = max(marker_size, 5)
            marker_size = min(marker_size, 50)
            self.settings.marker_size = marker_size
        self.settings_panel.marker_size_ctrl.SetValue(self.settings.marker_size)
        self.settings_panel.add_experiments_buttons()

    def OnActive(self, event):
        if self.IsShown() and type(self.viewer).__name__ != "_wxPyDeadObject":
            self.viewer.Refresh()

    def OnClose(self, event):
        self.Unbind(wx.EVT_ACTIVATE)
        self.Destroy()
        event.Skip()

    def OnDestroy(self, event):
        if self.parent is not None:
            self.parent.viewer = None
        event.Skip()

    def OnKeyDown(self, event):
        key = event.GetUnicodeKey()
        if key == wx.WXK_NONE:
            key = event.GetKeyCode()
        dxs = {wx.WXK_LEFT: -1, wx.WXK_RIGHT: +1, wx.WXK_UP: 0, wx.WXK_DOWN: 0}
        dys = {wx.WXK_LEFT: 0, wx.WXK_RIGHT: 0, wx.WXK_UP: +1, wx.WXK_DOWN: -1}

        if key in dxs:
            dx = dxs[key]
            dy = dys[key]
            if event.ShiftDown():
                scale = 0.1
            else:
                scale = 1.0
            self.do_Step(dx, dy, scale)

    def do_Step(self, dx, dy, scale):
        v = self.viewer
        rc = v.rotation_center
        gl.glMatrixMode(gl.GL_MODELVIEW)
        gltbx.util.rotate_object_about_eye_x_and_y(
            scale, rc[0], rc[1], rc[2], dx, dy, 0, 0
        )
        v.OnRedraw()

    def create_viewer_panel(self):
        if self.settings.black_background:
            background_rgb = (0, 0, 0)
        else:
            background_rgb = (255, 255, 255)
        self.viewer = RLVWindow(
            settings=self.settings,
            parent=self,
            size=(800, 600),
            background_rgb=background_rgb,
        )

    def create_settings_panel(self):
        self.settings_panel = SettingsWindow(self, -1, style=wx.RAISED_BORDER)

    def set_points(self):
        Render3d.set_points(self)
        self.settings_panel.d_min_ctrl.SetValue(self.settings.d_min)
        self.settings_panel.z_min_ctrl.SetValue(self.settings.z_min)
        self.settings_panel.z_max_ctrl.SetValue(self.settings.z_max)
        self.settings_panel.n_min_ctrl.SetValue(self.settings.n_min)
        self.settings_panel.n_max_ctrl.SetValue(self.settings.n_max)
        if self.settings.partiality_min is not None:
            self.settings_panel.partiality_min_ctrl.SetValue(
                self.settings.partiality_min
            )
        if self.settings.partiality_max is not None:
            self.settings_panel.partiality_max_ctrl.SetValue(
                self.settings.partiality_max
            )

    def update_settings(self, *args, **kwds):
        if self.settings.beam_centre:
            self.set_beam_centre(
                self.settings.beam_centre_panel, self.settings.beam_centre
            )
        self.map_points_to_reciprocal_space()
        self.set_points()
        self.viewer.update_settings(*args, **kwds)

    def update_statusbar(self):
        model_view_matrix = gltbx.util.get_gl_modelview_matrix()
        txt = (
            "Model view matrix: "
            + "["
            + ", ".join("%.4f" % m for m in model_view_matrix)
            + "]"
        )
        self.statusbar.SetStatusText(txt)


class SettingsWindow(wxtbx.utils.SettingsPanel):
    def __init__(self, *args, **kwds):
        wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
        self.Bind(wx.EVT_CHAR, self.OnChar)

    def OnChar(self, event):
        self.GetParent().viewer.OnChar(event)

    def add_controls(self):
        # d_min control

        self.d_min_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=0.05, min_val=0, digits=2
        )
        self.d_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.d_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "High resolution:")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.d_min_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.d_min_ctrl)

        self.z_min_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=1, min_val=0, digits=0
        )
        self.z_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.z_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Min Z")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.z_min_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.z_min_ctrl)

        self.z_max_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=1, min_val=0, digits=0
        )
        self.z_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.z_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Max Z")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.z_max_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.z_max_ctrl)

        # Control for spot size (utility depends on n_signal column in reflection
        # file - will be ignored if not in file

        self.n_min_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=1, min_val=0, digits=0
        )
        self.n_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.n_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Min Pixels")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.n_min_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.n_min_ctrl)

        self.n_max_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=1, min_val=0, digits=0
        )
        self.n_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.n_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Max Pixels")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.n_max_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.n_max_ctrl)

        # end new control

        self.partiality_min_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=0.01, digits=3, min_val=0, max_val=1
        )
        self.partiality_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.partiality_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Min partiality")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.partiality_min_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(
            floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.partiality_min_ctrl
        )

        self.partiality_max_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=0.01, digits=3, min_val=0, max_val=1
        )
        self.partiality_max_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.partiality_max_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Max partiality")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.partiality_max_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(
            floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.partiality_max_ctrl
        )

        ctrls = self.create_controls(
            setting="show_rotation_axis", label="Show rotation axis"
        )
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
        ctrls = self.create_controls(
            setting="show_beam_vector", label="Show beam vector"
        )
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
        ctrls = self.create_controls(
            setting="show_reciprocal_cell", label="Show reciprocal cell"
        )
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
        ctrls = self.create_controls(
            setting="label_nearest_point", label="Label nearest point"
        )
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
        self.reverse_phi_ctrl = self.create_controls(
            setting="reverse_phi", label="Invert rotation axis"
        )[0]
        self.panel_sizer.Add(self.reverse_phi_ctrl, 0, wx.ALL, 5)

        self.Bind(wx.EVT_CHECKBOX, self.OnChangeSettings, self.reverse_phi_ctrl)

        self.crystal_frame_tooltip = wx.ToolTip(
            "Show the reciprocal lattice(s) in the crystal rather than the laboratory frame"
        )
        self.crystal_frame_ctrl = self.create_controls(
            setting="crystal_frame", label="Show in crystal frame"
        )[0]
        self.crystal_frame_ctrl.SetToolTip(self.crystal_frame_tooltip)
        self.panel_sizer.Add(self.crystal_frame_ctrl, 0, wx.ALL, 5)

        self.Bind(wx.EVT_CHECKBOX, self.OnChangeSettings, self.crystal_frame_ctrl)

        self.beam_panel_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), min_val=0, increment=1, digits=0
        )
        self.beam_panel_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.beam_panel_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Beam centre panel")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.beam_panel_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.beam_panel_ctrl)

        self.beam_fast_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=0.01, digits=2
        )
        self.beam_fast_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.beam_fast_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        # box = wx.BoxSizer(wx.HORIZONTAL)
        grid = wx.FlexGridSizer(cols=2, rows=2, vgap=0, hgap=0)
        self.panel_sizer.Add(grid)
        label = wx.StaticText(self, -1, "Beam centre (mm)")
        grid.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        grid.Add(self.beam_fast_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.beam_fast_ctrl)

        self.beam_slow_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=0.01, digits=2
        )
        self.beam_slow_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.beam_slow_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        grid.Add(wx.StaticText(self, -1, ""), 0, wx.EXPAND)
        grid.Add(self.beam_slow_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.beam_slow_ctrl)

        self.marker_size_ctrl = floatspin.FloatSpin(
            parent=self, size=(150, -1), increment=1, digits=0, min_val=1
        )
        self.marker_size_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.marker_size_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Marker size:")
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(self.marker_size_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.marker_size_ctrl)

        self.btn = SegmentedRadioControl(self, style=SEGBTN_HORIZONTAL)
        self.btn.AddSegment("all")
        self.btn.AddSegment("indexed")
        self.btn.AddSegment("unindexed")
        self.btn.AddSegment("integrated")
        self.btn.SetSelection(
            ["all", "indexed", "unindexed", "integrated"].index(self.settings.display)
        )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnChangeSettings, self.btn)
        self.GetSizer().Add(self.btn, 0, wx.ALL, 5)

        self.outlier_btn = SegmentedRadioControl(self, style=SEGBTN_HORIZONTAL)
        self.outlier_btn.AddSegment("all")
        self.outlier_btn.AddSegment("inliers")
        self.outlier_btn.AddSegment("outliers")
        self.outlier_btn.SetSelection(
            [None, "inliers", "outliers"].index(self.settings.outlier_display)
        )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnChangeSettings, self.outlier_btn)
        self.GetSizer().Add(self.outlier_btn, 0, wx.ALL, 5)

    def add_value_widgets(self, sizer):
        sizer.Add(
            wx.StaticText(self.panel, -1, "Value:"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.value_info = wx.TextCtrl(
            self.panel, -1, size=(80, -1), style=wx.TE_READONLY
        )
        sizer.Add(self.value_info, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)

    def add_experiments_buttons(self):
        n = flex.max(self.parent.reflections_input["id"])
        if n <= 0 or n > self.settings.max_experiments:
            self.expt_btn = None
            return

        box = wx.BoxSizer(wx.VERTICAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Experiment ids:")
        box.Add(label, 0, wx.ALL, 5)

        self.expt_btn = SegmentedToggleControl(self, style=SEGBTN_HORIZONTAL)
        for i in range(-1, n + 1):
            self.expt_btn.AddSegment(str(i))
            if (
                self.settings.experiment_ids is not None
                and i in self.settings.experiment_ids
            ):
                self.expt_btn.SetValue(i + 1, True)

        self.expt_btn.Realize()
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnChangeSettings, self.expt_btn)
        box.Add(self.expt_btn, 0, wx.ALL, 5)

    def OnChangeSettings(self, event):
        self.settings.d_min = self.d_min_ctrl.GetValue()
        self.settings.z_min = self.z_min_ctrl.GetValue()
        self.settings.z_max = self.z_max_ctrl.GetValue()
        self.settings.n_min = int(self.n_min_ctrl.GetValue())
        self.settings.n_max = int(self.n_max_ctrl.GetValue())
        self.settings.partiality_min = self.partiality_min_ctrl.GetValue()
        self.settings.partiality_max = self.partiality_max_ctrl.GetValue()

        old_beam_panel = self.settings.beam_centre_panel
        old_beam_centre = self.settings.beam_centre
        self.settings.beam_centre_panel = self.beam_panel_ctrl.GetValue()
        if self.beam_fast_ctrl.IsEnabled() and self.beam_slow_ctrl.IsEnabled():
            self.settings.beam_centre = (
                self.beam_fast_ctrl.GetValue(),
                self.beam_slow_ctrl.GetValue(),
            )
        else:
            self.settings.beam_centre = None
        self.settings.reverse_phi = self.reverse_phi_ctrl.GetValue()
        self.settings.crystal_frame = self.crystal_frame_ctrl.GetValue()
        self.settings.marker_size = self.marker_size_ctrl.GetValue()
        for i, display in enumerate(("all", "indexed", "unindexed", "integrated")):
            if self.btn.values[i]:
                self.settings.display = display
                break

        for i, display in enumerate(("all", "inliers", "outliers")):
            if self.outlier_btn.values[i]:
                self.settings.outlier_display = display
                break

        if self.expt_btn is not None:
            expt_ids = []
            for i in range(len(self.expt_btn.segments)):
                if self.expt_btn.GetValue(i):
                    expt_ids.append(i - 1)
            self.settings.experiment_ids = expt_ids

        try:
            self.parent.update_settings()
        except ValueError:  # Handle beam centre changes, which could fail
            self.settings.beam_centre_panel = old_beam_panel
            self.settings.beam_centre = old_beam_centre
            self.beam_panel_ctrl.SetValue(old_beam_panel)
            self.beam_fast_ctrl.SetValue(old_beam_centre[0])
            self.beam_slow_ctrl.SetValue(old_beam_centre[1])


class RLVWindow(wx_viewer.show_points_and_lines_mixin):
    def __init__(self, settings, *args, **kwds):
        super().__init__(*args, **kwds)
        self.settings = settings
        self.points = flex.vec3_double()
        self.colors = None
        self.palette = None
        self.rotation_axis = None
        self.beam_vector = None
        self.recip_latt_vectors = None
        self.recip_crystal_vectors = None
        self.flag_show_minimum_covering_sphere = False
        self.minimum_covering_sphere = None
        self.field_of_view_y = 0.001
        if self.settings.autospin:
            self.autospin_allowed = True
            self.yspin = 1
            self.xspin = 1
            self.autospin = True

    def set_points(self, points):
        self.points = points
        self.points_display_list = None
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()

    def set_points_data(self, reflections):
        dstar = reflections["rlp"].norms()
        dstar.set_selected(dstar == 0, 1e-8)
        self.points_data = {
            "panel": reflections["panel"],
            "id": reflections["id"],
            "xyz": reflections["xyzobs.px.value"],
            "d_spacing": 1 / dstar,
        }
        if "miller_index" in reflections:
            self.points_data["miller_index"] = reflections["miller_index"]

    def set_colors(self, colors):
        assert len(colors) == len(self.points)
        self.colors = colors

    def set_palette(self, palette):
        self.palette = palette

    def draw_points(self):
        if self.points_display_list is None:
            self.points_display_list = gltbx.gl_managed.display_list()
            self.points_display_list.compile()
            gl.glLineWidth(1)
            if self.colors is None:
                self.colors = flex.vec3_double(len(self.points), (1, 1, 1))
            for point, color in zip(self.points, self.colors):
                self.draw_cross_at(point, color=color)
            self.points_display_list.end()
        self.points_display_list.call()

    def set_rotation_axis(self, axis):
        self.rotation_axis = axis

    def set_beam_vector(self, beam):
        self.beam_vector = beam

    def set_reciprocal_lattice_vectors(self, vectors_per_crystal):
        self.recip_latt_vectors = vectors_per_crystal

    def set_reciprocal_crystal_vectors(self, vectors_per_crystal):
        self.recip_crystal_vectors = vectors_per_crystal

    # --- user input and settings
    def update_settings(self):
        self.points_display_list = None
        self.Refresh()

    def update_minimum_covering_sphere(self):
        n_points = min(1000, self.points.size())
        isel = flex.random_permutation(self.points.size())[:n_points]
        self.minimum_covering_sphere = minimum_covering_sphere(self.points.select(isel))

    def draw_cross_at(self, xyz, color=(1, 1, 1), f=None):
        (x, y, z) = xyz
        if f is None:
            f = 0.01 * self.settings.marker_size
        wx_viewer.show_points_and_lines_mixin.draw_cross_at(
            self, (x, y, z), color=color, f=f
        )

    def DrawGL(self):
        wx_viewer.show_points_and_lines_mixin.DrawGL(self)
        if self.rotation_axis is not None and self.settings.show_rotation_axis:
            self.draw_axis(self.rotation_axis, "phi")
        if self.beam_vector is not None and self.settings.show_beam_vector:
            self.draw_axis(self.beam_vector, "beam")

        if self.settings.show_reciprocal_cell:
            # if we don't have one sort of vector we don't have the other either
            vectors = self.recip_latt_vectors
            if self.settings.crystal_frame:
                vectors = self.recip_crystal_vectors

            if vectors:
                for i, axes in enumerate(vectors):
                    if self.settings.experiment_ids:
                        if i not in self.settings.experiment_ids:
                            continue
                    j = (i + 1) % self.palette.size()
                    color = self.palette[j]
                    self.draw_cell(axes, color)

        if self.settings.label_nearest_point:
            self.label_nearest_point()

        self.GetParent().update_statusbar()

    def draw_axis(self, axis, label):
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()
        s = self.minimum_covering_sphere
        scale = max(max(s.box_max()), abs(min(s.box_min())))
        gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
        gl.glDisable(gl.GL_LIGHTING)
        if self.settings.black_background:
            gl.glColor3f(1.0, 1.0, 1.0)
        else:
            gl.glColor3f(0.0, 0.0, 0.0)
        gl.glLineWidth(1.0)
        gl.glBegin(gl.GL_LINES)
        gl.glVertex3f(0.0, 0.0, 0.0)
        gl.glVertex3f(axis[0] * scale, axis[1] * scale, axis[2] * scale)
        gl.glEnd()
        gl.glRasterPos3f(
            0.5 + axis[0] * scale, 0.2 + axis[1] * scale, 0.2 + axis[2] * scale
        )
        gltbx.fonts.ucs_bitmap_8x13.render_string(label)
        gl.glEnable(gl.GL_LINE_STIPPLE)
        gl.glLineStipple(4, 0xAAAA)
        gl.glBegin(gl.GL_LINES)
        gl.glVertex3f(0.0, 0.0, 0.0)
        gl.glVertex3f(-axis[0] * scale, -axis[1] * scale, -axis[2] * scale)
        gl.glEnd()
        gl.glDisable(gl.GL_LINE_STIPPLE)

    def draw_cell(self, axes, color):
        astar, bstar, cstar = axes[0], axes[1], axes[2]
        gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
        gl.glDisable(gl.GL_LIGHTING)
        gl.glColor3f(*color)
        gl.glLineWidth(2.0)
        gl.glBegin(gl.GL_LINES)
        gl.glVertex3f(0.0, 0.0, 0.0)
        gl.glVertex3f(*astar.elems)
        gl.glVertex3f(0.0, 0.0, 0.0)
        gl.glVertex3f(*bstar.elems)
        gl.glVertex3f(0.0, 0.0, 0.0)
        gl.glVertex3f(*cstar.elems)
        gl.glEnd()
        gl.glRasterPos3f(*(1.01 * astar).elems)
        gltbx.fonts.ucs_bitmap_8x13.render_string("a*")
        gl.glRasterPos3f(*(1.01 * bstar).elems)
        gltbx.fonts.ucs_bitmap_8x13.render_string("b*")
        gl.glRasterPos3f(*(1.01 * cstar).elems)
        gltbx.fonts.ucs_bitmap_8x13.render_string("c*")
        gl.glEnable(gl.GL_LINE_STIPPLE)
        gl.glLineStipple(4, 0xAAAA)
        farpoint = astar + bstar + cstar
        # a* face
        gl.glBegin(gl.GL_LINE_LOOP)
        gl.glVertex3f(*farpoint.elems)
        gl.glVertex3f(*(farpoint - bstar).elems)
        gl.glVertex3f(*(farpoint - bstar - cstar).elems)
        gl.glVertex3f(*(farpoint - cstar).elems)
        gl.glEnd()
        # b* face
        gl.glBegin(gl.GL_LINE_LOOP)
        gl.glVertex3f(*farpoint.elems)
        gl.glVertex3f(*(farpoint - astar).elems)
        gl.glVertex3f(*(farpoint - astar - cstar).elems)
        gl.glVertex3f(*(farpoint - cstar).elems)
        gl.glEnd()
        # c* face
        gl.glBegin(gl.GL_LINE_LOOP)
        gl.glVertex3f(*farpoint.elems)
        gl.glVertex3f(*(farpoint - bstar).elems)
        gl.glVertex3f(*(farpoint - bstar - astar).elems)
        gl.glVertex3f(*(farpoint - astar).elems)
        gl.glEnd()
        gl.glDisable(gl.GL_LINE_STIPPLE)

    def label_nearest_point(self):
        ann = AnnAdaptorSelfInclude(self.points.as_double(), 3)
        ann.query(self.rotation_center)
        i = ann.nn[0]
        gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
        gl.glDisable(gl.GL_LIGHTING)
        gl.glColor3f(1.0, 1.0, 1.0)
        gl.glLineWidth(1.0)
        xyz = self.points_data["xyz"][i]
        exp_id = self.points_data["id"][i]
        panel = self.points_data["panel"][i]
        d_spacing = self.points_data["d_spacing"][i]
        label = (
            f"id: {exp_id}; panel: {panel}\n"
            f"xyz: {xyz[0]:.1f} {xyz[1]:.1f} {xyz[2]:.1f}\n"
            f"res: {d_spacing:.2f} Angstrom"
        )
        if "miller_index" in self.points_data and exp_id != -1:
            hkl = self.points_data["miller_index"][i]
            label += f"\nhkl: {hkl}"
        line_spacing = round(gltbx.fonts.ucs_bitmap_8x13.height())
        for j, string in enumerate(label.splitlines()):
            gl.glRasterPos3f(*self.points[i])
            gl.glBitmap(0, 0, 0.0, 0.0, line_spacing, -j * line_spacing, b" ")
            gltbx.fonts.ucs_bitmap_8x13.render_string(string)

    def rotate_view(self, x1, y1, x2, y2, shift_down=False, scale=0.1):
        super().rotate_view(x1, y1, x2, y2, shift_down=shift_down, scale=scale)

    def OnLeftUp(self, event):
        self.was_dragged = True
        super().OnLeftUp(event)

    def initialize_modelview(self, eye_vector=None, angle=None):
        super().initialize_modelview(eye_vector=eye_vector, angle=angle)
        self.rotation_center = (0, 0, 0)
        self.move_to_center_of_viewport(self.rotation_center)
        if self.settings.model_view_matrix is not None:
            gl.glLoadMatrixd(self.settings.model_view_matrix)
