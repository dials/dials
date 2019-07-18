# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import absolute_import, division, print_function
from dials.util import wx_viewer
import copy
import wx
import wxtbx.utils
from gltbx.gl import *
import gltbx
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
import libtbx.phil
from dxtbx.model import MultiAxisGoniometer

help_message = """

"""

phil_scope = libtbx.phil.parse(
    """
  angle = None
    .type = float
    .multiple = True
  detector_distance = None
    .type = float
  z_min = None
    .type = float
  z_max = None
    .type = float
  marker_size = 100
    .type = int(value_min=1)
  autospin = False
    .type = bool
  predict = False
    .type = bool
  prediction_width = None
    .type = float(value_min=0)
    .help = "Width of prediction window (degrees)"
  show_panel_axes = False
    .type = bool
    .help = "Plot the fast, slow and normal vectors for each panel."
  show_crystal_axes = True
    .type = bool
    .help = "Plot the crystal reciprocal cell axes"
  require_images = True
    .type = bool
    .help = "Flag which can be set to False to launch image viewer without
      checking the image format (needed for some image format classes).
      Alternative to DIALS_EXPORT_DO_NOT_CHECK_FORMAT environment variable."
"""
)


def settings():
    return phil_scope.fetch().extract()


class render_3d(object):
    def __init__(self):
        self.goniometer_orig = None
        self.crystal = None

    def load_imageset(self, imageset, crystal=None):
        self.imageset = imageset
        if crystal is not None:
            self.crystal = crystal
        if self.imageset.get_goniometer() is not None:
            from scitbx import matrix

            gonio = self.imageset.get_goniometer()
            axis = matrix.col(gonio.get_rotation_axis()).elems
            self.viewer.set_rotation_axis(axis)
        self.viewer.set_beam_vector(self.imageset.get_beam().get_s0())

        self.set_points()

    def set_goniometer_points(self):
        gonio = self.imageset.get_goniometer()
        detector = self.imageset.get_detector()
        beam = self.imageset.get_beam()

        angle = self.settings.angle
        if angle:
            assert len(angle) == len(gonio.get_angles())
            gonio.set_angles(angle)

        distance = self.settings.detector_distance
        if distance:
            import math
            from scitbx import matrix

            p_id = detector.get_panel_intersection(beam.get_s0())
            if p_id >= 0:
                if len(detector) > 1:
                    p = detector.hierarchy()
                else:
                    p = detector[0]
                d_normal = matrix.col(p.get_normal())
                d_origin = matrix.col(p.get_origin())
                d_distance = math.fabs(
                    d_origin.dot(d_normal) - p.get_directed_distance()
                )
                assert d_distance < 0.001, d_distance
                translation = d_normal * (distance - p.get_directed_distance())
                new_origin = d_origin + translation
                d_distance = math.fabs(new_origin.dot(d_normal) - distance)
                assert d_distance < 0.001, d_distance
                p.set_frame(p.get_fast_axis(), p.get_slow_axis(), new_origin.elems)

        gonio_masker = (
            self.imageset.get_format_class()
            .get_instance(self.imageset.paths()[0], **self.imageset.data().get_params())
            .get_masker(goniometer=gonio)
        )
        if gonio_masker is None:
            return

        points = gonio_masker.extrema_at_scan_angle(
            gonio.get_angles()[gonio.get_scan_axis()]
        )
        points.insert(0, (0, 0, 0))

        line_i_seqs = flex.vec2_double((0, i) for i in range(1, points.size()))
        line_i_seqs += (self.viewer.points.size(), self.viewer.points.size())
        for i_seqs in line_i_seqs:
            self.viewer.line_i_seqs.append([int(i_seq) for i_seq in i_seqs])
            self.viewer.line_colors[(i_seqs)] = (100 / 255, 120 / 255, 255 / 255)

        self.viewer.points.extend(points)

        shadow = gonio_masker.project_extrema(
            detector, gonio.get_angles()[gonio.get_scan_axis()]
        )

        for shadow_points, p in zip(shadow, detector):
            n = self.viewer.points.size()
            line_i_seqs = []
            line_colors = {}
            for i in range(shadow_points.size()):
                if i < shadow_points.size() - 1:
                    line_i_seqs.append((n + i, n + i + 1))
                else:
                    line_i_seqs.append((n + i, n))
                line_colors[line_i_seqs[-1]] = (1, 1, 1)
            self.viewer.line_colors.update(line_colors)
            self.viewer.line_i_seqs.extend(line_i_seqs)
            self.viewer.points.extend(
                p.get_lab_coord(shadow_points * p.get_pixel_size()[0])
            )

    def set_detector_points(self):
        detector = self.imageset.get_detector()
        points = flex.vec3_double()
        line_i_seqs = flex.vec2_double()
        i = 0
        for p in detector:
            image_size = p.get_image_size_mm()
            points.append(p.get_lab_coord((0, 0)))
            points.append(p.get_lab_coord((0, image_size[1])))
            points.append(p.get_lab_coord(image_size))
            points.append(p.get_lab_coord((image_size[0], 0)))
            line_i_seqs.append((i, i + 1))
            line_i_seqs.append((i + 1, i + 2))
            line_i_seqs.append((i + 2, i + 3))
            line_i_seqs.append((i + 3, i))
            i += 4
        line_i_seqs += (self.viewer.points.size(), self.viewer.points.size())
        self.viewer.points.extend(points)
        for i_seqs in line_i_seqs:
            self.viewer.line_i_seqs.append([int(i_seq) for i_seq in i_seqs])

    def set_points(self):
        # reset point/line lists
        self.viewer.reset()

        try:
            self.set_goniometer_points()
        except Exception:
            pass
        self.set_detector_points()
        self.viewer.update_minimum_covering_sphere()
        if self.settings.predict and self.crystal is not None:
            self.set_reflection_points()

    def set_reflection_points(self):
        import time

        t0 = time.time()
        predicted = self.predict()
        if predicted is None:
            return
        t1 = time.time()
        print("Predicted %i reflections in %.2f s" % (predicted.size(), (t1 - t0)))
        xc, yc, zc = predicted["xyzcal.mm"].parts()
        xycal = flex.vec2_double(xc, yc)
        panel = predicted["panel"]
        detector = self.imageset.get_detector()
        for pid, p in enumerate(detector):
            self.viewer.points.extend(p.get_lab_coord(xycal.select(panel == pid)))

    def predict(self):
        assert self.crystal is not None
        from dxtbx.model.experiment_list import Experiment

        imageset = self.imageset
        scan = copy.deepcopy(imageset.get_scan())
        gonio = imageset.get_goniometer()
        prediction_width = self.settings.prediction_width
        if prediction_width is None:
            prediction_width = scan.get_oscillation()[1]
        if isinstance(gonio, MultiAxisGoniometer):
            scan_angle = gonio.get_angles()[gonio.get_scan_axis()]
        else:
            return
        scan.set_oscillation((scan_angle, prediction_width))
        expt = Experiment(
            imageset=imageset,
            crystal=self.crystal,
            detector=imageset.get_detector(),
            beam=imageset.get_beam(),
            scan=scan[:1],
            goniometer=imageset.get_goniometer(),
        )

        # Populate the reflection table with predictions
        from dials.array_family import flex

        predicted = flex.reflection_table.from_predictions(expt, force_static=True)
        predicted["id"] = flex.int(len(predicted), 0)
        return predicted


class ExperimentViewer(wx.Frame, render_3d):
    def __init__(self, *args, **kwds):
        wx.Frame.__init__(self, *args, **kwds)
        render_3d.__init__(self)
        self.parent = self.GetParent()
        self.statusbar = self.CreateStatusBar()
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        app = wx.GetApp()
        if getattr(app, "settings", None) is not None:
            # XXX copying the initial settings avoids awkward interactions when
            # multiple viewer windows are opened
            self.settings = copy.deepcopy(app.settings)
        else:
            self.settings = settings()

        self.create_settings_panel()
        self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
        self.create_viewer_panel()
        self.sizer.Add(self.viewer, 1, wx.EXPAND | wx.ALL)
        self.SetSizer(self.sizer)
        self.sizer.SetSizeHints(self)
        self.Bind(wx.EVT_CLOSE, self.OnClose, self)
        self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
        self.Bind(wx.EVT_ACTIVATE, self.OnActive)
        self.viewer.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
        self.viewer.SetFocus()

    def load_imageset(self, imageset, crystal=None):
        render_3d.load_imageset(self, imageset, crystal=crystal)
        self.settings_panel.add_goniometer_controls(imageset.get_goniometer())

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
        glMatrixMode(GL_MODELVIEW)
        gltbx.util.rotate_object_about_eye_x_and_y(
            scale, rc[0], rc[1], rc[2], dx, dy, 0, 0
        )
        v.OnRedraw()

    def create_viewer_panel(self):
        self.viewer = GeometryWindow(
            settings=self.settings,
            parent=self,
            size=(800, 600),
            # orthographic=True
        )

    def create_settings_panel(self):
        self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

    def set_points(self):
        render_3d.set_points(self)

    def update_settings(self, *args, **kwds):
        self.set_points()
        self.viewer.update_settings(*args, **kwds)


class settings_window(wxtbx.utils.SettingsPanel):
    def __init__(self, *args, **kwds):
        wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
        self.Bind(wx.EVT_CHAR, self.OnChar)

    def OnChar(self, event):
        self.GetParent().viewer.OnChar(event)

    def add_controls(self):

        ctrls = self.create_controls(setting="show_panel_axes", label="Show panel axes")
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)

        ctrls = self.create_controls(
            setting="show_crystal_axes", label="Show crystal axes"
        )
        self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)

    def add_goniometer_controls(self, goniometer):
        from wx.lib.agw import floatspin

        self.distance_ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=2)
        self.distance_ctrl.SetValue(self.settings.detector_distance)
        self.distance_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
        if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
            self.distance_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_sizer.Add(box)
        label = wx.StaticText(self, -1, "Detector distance")
        box.Add(self.distance_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, self.distance_ctrl)

        self.rotation_angle_ctrls = []
        if isinstance(goniometer, MultiAxisGoniometer):
            names = goniometer.get_names()
            axes = goniometer.get_axes()
            angles = goniometer.get_angles()
            for name, axis, angle in zip(names, axes, angles):
                ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=3)
                ctrl.SetValue(angle)
                ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
                if wx.VERSION >= (2, 9):  # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
                    ctrl.SetBackgroundColour(self.GetBackgroundColour())
                box = wx.BoxSizer(wx.HORIZONTAL)
                self.panel_sizer.Add(box)
                label = wx.StaticText(self, -1, "%s angle" % name)
                box.Add(ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
                box.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
                self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, ctrl)
                self.rotation_angle_ctrls.append(ctrl)

    def OnChangeSettings(self, event):
        self.settings.detector_distance = self.distance_ctrl.GetValue()
        for i, ctrl in enumerate(self.rotation_angle_ctrls):
            self.settings.angle[i] = ctrl.GetValue()
        self.parent.update_settings()


class GeometryWindow(wx_viewer.show_points_and_lines_mixin):
    def __init__(self, settings, *args, **kwds):
        super(GeometryWindow, self).__init__(*args, **kwds)
        self.settings = settings
        self.points = flex.vec3_double()
        self.colors = None
        self.rotation_axis = None
        self.beam_vector = None
        self.flag_show_minimum_covering_sphere = False
        self.minimum_covering_sphere = None
        self.field_of_view_y = 20
        if self.settings.autospin:
            self.autospin_allowed = True
            self.yspin = 1
            self.xspin = 1
            self.autospin = True

    def reset(self):
        self.labels = []
        self.points = flex.vec3_double()
        self.line_i_seqs = []
        self.line_colors = {}
        self.spheres = []
        self.labels_display_list = None
        self.points_display_list = None
        self.lines_display_list = None

    def set_points(self, points):
        self.points = points
        self.points_display_list = None
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()

    def set_colors(self, colors):
        assert len(colors) == len(self.points)
        self.colors = colors

    def draw_points(self):
        if self.points_display_list is None:
            self.points_display_list = gltbx.gl_managed.display_list()
            self.points_display_list.compile()
            glLineWidth(1)
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

    # --- user input and settings
    def update_settings(self):
        self.points_display_list = None
        self.Refresh()

    def update_minimum_covering_sphere(self):
        self.minimum_covering_sphere = minimum_covering_sphere(
            self.points, epsilon=1e-3
        )

    def draw_cross_at(self, xyz, color=(1, 1, 1), f=None):
        (x, y, z) = xyz
        if f is None:
            f = 0.01 * self.settings.marker_size
        wx_viewer.show_points_and_lines_mixin.draw_cross_at(
            self, (x, y, z), color=color, f=f
        )

    def DrawGL(self):
        wx_viewer.show_points_and_lines_mixin.DrawGL(self)
        gonio = self.parent.imageset.get_goniometer()
        beam = self.parent.imageset.get_beam()
        detector = self.parent.imageset.get_detector()
        if self.settings.show_panel_axes and detector is not None:
            for p in detector:
                image_size = p.get_image_size_mm()
                from scitbx import matrix

                o = matrix.col(p.get_lab_coord((0, 0)))
                f = matrix.col(p.get_lab_coord((image_size[0], 0)))
                _f = (f - o).normalize()
                s = matrix.col(p.get_lab_coord((0, image_size[1])))
                _s = (s - o).normalize()
                _n = _f.cross(_s)
                n = _n * (0.25 * (f.length() + s.length())) + o
                self.draw_lab_axis(o.elems, f.elems, "FAST")
                self.draw_lab_axis(o.elems, s.elems, "SLOW")
                self.draw_lab_axis(o.elems, n.elems, "NORM")

        from scitbx import matrix

        if isinstance(gonio, MultiAxisGoniometer):
            R = matrix.identity(3)
            names = reversed(gonio.get_names())
            axes = reversed(gonio.get_axes())
            angles = reversed(gonio.get_angles())
            for name, axis, angle in zip(names, axes, angles):
                axis = R * matrix.col(axis)
                self.draw_axis(axis.elems, name)
                R = axis.axis_and_angle_as_r3_rotation_matrix(angle, deg=True) * R
        elif gonio is not None:
            axis = matrix.col(gonio.get_rotation_axis())
            self.draw_axis(axis.elems, "phi")
        self.draw_axis(beam.get_s0(), "beam")
        crystal = self.parent.crystal
        if self.settings.show_crystal_axes and crystal is not None:
            crystal = copy.deepcopy(crystal)
            scan = self.parent.imageset.get_scan()
            fixed_rotation = matrix.sqr(gonio.get_fixed_rotation())
            setting_rotation = matrix.sqr(gonio.get_setting_rotation())
            rotation_axis = matrix.col(gonio.get_rotation_axis_datum())
            if isinstance(gonio, MultiAxisGoniometer):
                angle = gonio.get_angles()[gonio.get_scan_axis()]
            else:
                angle = scan.get_oscillation()[0]
            rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
                angle, deg=True
            )
            U0 = matrix.sqr(crystal.get_U())

            # Goniometer datum setting [D] at which the orientation was determined
            D = (setting_rotation * rotation_matrix * fixed_rotation).inverse()

            U = D.inverse() * U0
            B = matrix.sqr(crystal.get_B())
            a_star = U * B * matrix.col((1, 0, 0))
            b_star = U * B * matrix.col((0, 1, 0))
            c_star = U * B * matrix.col((0, 0, 1))
            color = (1.0, 0.0, 0.0)  # red
            self.draw_axis(a_star.normalize() * 0.5, "a*", color=color)
            self.draw_axis(b_star.normalize() * 0.5, "b*", color=color)
            self.draw_axis(c_star.normalize() * 0.5, "c*", color=color)

    def draw_axis(self, axis, label, color=(1.0, 1.0, 1.0)):
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()
        s = self.minimum_covering_sphere
        scale = max(max(s.box_max()), abs(min(s.box_min())))
        gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glLineWidth(1.0)
        glBegin(GL_LINES)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(axis[0] * scale, axis[1] * scale, axis[2] * scale)
        glEnd()
        glRasterPos3f(
            0.5 + axis[0] * scale, 0.2 + axis[1] * scale, 0.2 + axis[2] * scale
        )
        gltbx.fonts.ucs_bitmap_8x13.render_string(label)
        glEnable(GL_LINE_STIPPLE)
        glLineStipple(4, 0xAAAA)
        glBegin(GL_LINES)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(-axis[0] * scale, -axis[1] * scale, -axis[2] * scale)
        glEnd()
        glDisable(GL_LINE_STIPPLE)

    def draw_lab_axis(self, start, end, label):
        mid = tuple([0.5 * (s + e) for s, e in zip(start, end)])
        gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 0.0)
        glLineWidth(1.0)
        glBegin(GL_LINES)
        glVertex3f(*start)
        glVertex3f(*end)
        glEnd()
        glRasterPos3f(*mid)
        gltbx.fonts.ucs_bitmap_8x13.render_string(label)

    def rotate_view(self, x1, y1, x2, y2, shift_down=False, scale=0.1):
        super(GeometryWindow, self).rotate_view(
            x1, y1, x2, y2, shift_down=shift_down, scale=scale
        )

    def OnLeftUp(self, event):
        self.was_dragged = True
        super(GeometryWindow, self).OnLeftUp(event)

    def initialize_modelview(self, eye_vector=None, angle=None):
        super(GeometryWindow, self).initialize_modelview(
            eye_vector=eye_vector, angle=angle
        )
        self.rotation_center = (0, 0, 0)
        self.move_to_center_of_viewport(self.rotation_center)


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    import os

    usage = "dials.geometry_viewer [options] models.expt"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(quick_parse=True, show_diff_phil=True)

    if "DIALS_EXPORT_DO_NOT_CHECK_FORMAT" in os.environ:
        print(
            "\nWarning: use of DIALS_EXPORT_DO_NOT_CHECK_FORMAT environment variable"
            "\nis no longer recommended, the recommended command is require_images=False"
        )
    if not params.require_images or "DIALS_EXPORT_DO_NOT_CHECK_FORMAT" in os.environ:
        check_format = False
    else:
        check_format = True

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=check_format,
        epilog=help_message,
    )

    try:
        params, options = parser.parse_args(show_diff_phil=False)
    except Exception as e:
        print(e)
        print(
            "dials.geometry_viewer help: Error in parsing data may be due to missing \n"
            "files. If so, this may be overcome by require_images=False\n"
        )
        sys.exit()

    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    imagesets = experiments.imagesets()
    crystal = experiments[0].crystal

    assert len(imagesets) == 1, len(imagesets)
    imageset = imagesets[0]
    gonio = imageset.get_goniometer()
    if not params.detector_distance:
        detector = imageset.get_detector()
        if len(detector) > 1:
            params.detector_distance = detector.hierarchy().get_directed_distance()
        else:
            params.detector_distance = detector[0].get_directed_distance()

    if gonio is not None and not isinstance(gonio, MultiAxisGoniometer):
        from dxtbx.model.goniometer import GoniometerFactory

        gonio = GoniometerFactory.multi_axis(
            axes=flex.vec3_double((gonio.get_rotation_axis(),)),
            angles=flex.double((0,)),
            names=flex.std_string(("GON_OMEGA",)),
            scan_axis=0,
        )
        imageset.set_goniometer(gonio)

    if isinstance(gonio, MultiAxisGoniometer):
        if params.angle:
            assert len(params.angle) == len(gonio.get_angles())
        else:
            for angle in gonio.get_angles():
                params.angle.append(angle)

    import wxtbx.app

    a = wxtbx.app.CCTBXApp(0)
    a.settings = params
    f = ExperimentViewer(None, -1, "Experiment viewer", size=(1024, 768))
    f.load_imageset(imageset, crystal=crystal)
    f.Show()
    a.SetTopWindow(f)
    a.MainLoop()


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
