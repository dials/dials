from __future__ import annotations

import math
import json
import copy

import wx
from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

import cctbx.miller
from cctbx.crystal import symmetry
from scitbx.matrix import col


class UCSettingsFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        szr = wx.BoxSizer(wx.VERTICAL)
        self.phil_params = args[0].params
        panel = UCSettingsPanel(self)
        self.SetSizer(szr)
        szr.Add(panel, 1, wx.EXPAND)
        szr.Fit(panel)
        self.panel = panel
        self.sizer = szr
        self.Fit()
        self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)


class UCSettingsPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

        self.phil_params = args[0].phil_params

        # Needed to draw and delete the rings.  XXX Applies to
        # calibration_frame as well?
        self._pyslip = self.GetParent().GetParent().pyslip

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        # Number of decimal digits for distances.
        self.digits = 2

        # Wavelength control.
        beam = self._pyslip.tiles.raw_image.get_beam()
        self._wavelength = beam.get_wavelength()

        # Unit cell controls.
        if self.phil_params.calibrate_unit_cell.unit_cell is not None:
            self._cell = list(
                self.phil_params.calibrate_unit_cell.unit_cell.parameters()
            )
        else:
            self._cell = [4.18, 4.72, 58.38, 89.44, 89.63, 75.85]

        if self.phil_params.calibrate_unit_cell.space_group is not None:
            self._spacegroup = self.phil_params.calibrate_unit_cell.space_group
        else:
            self._spacegroup = "P1"

        self._show_hkl = self.phil_params.calibrate_unit_cell.show_hkl

        self._cell_control_names = [
            "uc_a_ctrl",
            "uc_b_ctrl",
            "uc_c_ctrl",
            "uc_alpha_ctrl",
            "uc_beta_ctrl",
            "uc_gamma_ctrl",
        ]

        box = wx.BoxSizer(wx.HORIZONTAL)

        self.uc_a = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[0],
            value=self._cell[0],
        )
        box.Add(
            self.uc_a, 0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5
        )
        box.Add(wx.StaticText(self, label="a"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_a)

        self.uc_alpha = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[3],
            value=self._cell[3],
        )
        box.Add(
            self.uc_alpha,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="alpha"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_alpha)

        sizer.Add(box)

        box = wx.BoxSizer(wx.HORIZONTAL)

        self.uc_b = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[1],
            value=self._cell[1],
        )
        box.Add(
            self.uc_b, 0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5
        )
        box.Add(wx.StaticText(self, label="b"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_b)

        self.uc_beta = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[4],
            value=self._cell[4],
        )
        box.Add(
            self.uc_beta, 0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5
        )
        box.Add(
            wx.StaticText(self, label="beta"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_beta)

        sizer.Add(box)

        box = wx.BoxSizer(wx.HORIZONTAL)

        self.uc_c = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[2],
            value=self._cell[2],
        )
        box.Add(
            self.uc_c, 0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5
        )
        box.Add(wx.StaticText(self, label="c"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_c)

        self.uc_gamma = FloatSpin(
            self,
            digits=self.digits,
            name=self._cell_control_names[5],
            value=self._cell[5],
        )
        box.Add(
            self.uc_gamma,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="gamma"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_gamma)

        sizer.Add(box)

        # Space group control
        box = wx.BoxSizer(wx.HORIZONTAL)

        self.space_group_ctrl = wx.TextCtrl(
            self, name="space group", value=self._spacegroup
        )
        box.Add(
            self.space_group_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="Space group"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(wx.EVT_TEXT, self.OnSpaceGroup, self.space_group_ctrl)

        sizer.Add(box)

        # Distance control
        img = self.GetParent().GetParent()._img
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.distance_ctrl = FloatSpin(
            self,
            digits=self.digits,
            name="Detector Distance",
            value=0,
        )
        self.distance_ctrl.SetIncrement(0.5)
        box.Add(
            self.distance_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )

        txtd = wx.StaticText(self, label="Detector Distance")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpin, self.distance_ctrl)
        sizer.Add(box)

        # Wavelength control
        img = self.GetParent().GetParent()._img
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.wavelength_ctrl = FloatSpin(
            self, digits=4, name="Wavelength", value=img.get_wavelength()
        )
        self.wavelength_ctrl.SetIncrement(0.05)
        box.Add(
            self.wavelength_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )

        txtw = wx.StaticText(self, label="Wavelength")
        box.Add(txtw, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpin, self.wavelength_ctrl)
        sizer.Add(box)

        # d_min control
        if self.phil_params.calibrate_unit_cell.d_min is not None:
            self.d_min = self.phil_params.calibrate_unit_cell.d_min
        else:
            self.d_min = 3.5
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.d_min_ctrl = FloatSpin(
            self, digits=self.digits, name="d_min", value=self.d_min
        )
        box.Add(
            self.d_min_ctrl,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )

        txtd = wx.StaticText(self, label="Highest resolution for ring display")
        box.Add(txtd, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpin, self.d_min_ctrl)
        sizer.Add(box)

        # Centering controls.
        self._center = [0, 0]
        box = wx.BoxSizer(wx.HORIZONTAL)

        self.spinner_fast = FloatSpin(
            self, digits=self.digits, name="fast_ctrl", value=self._center[0]
        )
        box.Add(
            self.spinner_fast,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="Center fast"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinCenter, self.spinner_fast)

        self.spinner_slow = FloatSpin(
            self, digits=self.digits, name="slow_ctrl", value=self._center[1]
        )
        box.Add(
            self.spinner_slow,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="Center slow"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        sizer.Add(box)

        # Rotation controls
        self._axis_rots = [0., 0.]
        box = wx.BoxSizer(wx.HORIZONTAL)

        self.spinner_rot_fast = FloatSpin(
            self, digits=self.digits, name="rot_fast_ctrl", value=self._axis_rots[0]
        )
        box.Add(
            self.spinner_rot_fast,
            0,
            wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="Rot fast"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinRot, self.spinner_rot_fast)

        self.spinner_rot_slow = FloatSpin(
            self, digits=self.digits, name="rot_slow_ctrl", value=self._axis_rots[1]
        )
        box.Add(
            self.spinner_rot_slow,
            0,
            wx.LEFT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        box.Add(
            wx.StaticText(self, label="Rot slow"),
            0,
            wx.ALL | wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        self.Bind(EVT_FLOATSPIN, self.OnSpinRot, self.spinner_rot_slow)
        sizer.Add(box)

        box = wx.BoxSizer(wx.HORIZONTAL)
        self.clear_button = wx.Button(self, -1, "Clear")
        box.Add(self.clear_button, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnClear, self.clear_button)
        sizer.Add(box)

        origin_box = wx.BoxSizer(wx.HORIZONTAL)
        self.origin = wx.StaticText(self, label="")
        origin_box.Add(self.origin, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(EVT_FLOATSPIN, self.OnSpinCenter, self.spinner_slow)

        detector = self._pyslip.tiles.raw_image.get_detector()
        d_biggest = 0
        for node in detector:
          d = node.get_origin()[2]
          if abs(d) > abs(d_biggest):
            d_biggest = d
        self.origin_cache = tuple(
            a+b for a,b in zip(detector.hierarchy().get_origin(), (0,0,d_biggest))
        )
        origin = detector.hierarchy().get_origin()
        new_origin = tuple(
            a+b for a,b in zip(origin, (0,0,d_biggest))
        )
        detector.hierarchy().set_frame(
            detector.hierarchy().get_fast_axis(),
            detector.hierarchy().get_slow_axis(),
            new_origin
        )

        for node in detector:
          origin = node.get_origin()
          new_origin = tuple(
              a-b for a,b in zip(origin, (0,0,d_biggest))
          )
          node.set_frame(
              node.get_fast_axis(),
              node.get_slow_axis(),
              new_origin
          )

        self.origin_cache = detector.hierarchy().get_origin()
        self.fast_cache = detector.hierarchy().get_fast_axis()
        self.slow_cache = detector.hierarchy().get_slow_axis()
        sizer.Add(origin_box)

        self.DrawRings()

    def __del__(self):
        if hasattr(self, "_ring_layer") and self._ring_layer is not None:
            self._pyslip.DeleteLayer(self._ring_layer)
            self._ring_layer = None

    def OnSpinCenter(self, event):
        obj = event.EventObject
        name = obj.GetName()

        if name == "fast_ctrl":
            self._center[0] = obj.GetValue()
        elif name == "slow_ctrl":
            self._center[1] = obj.GetValue()

        self.DrawRings()

    def OnSpinRot(self, event):
      obj = event.EventObject
      name = obj.GetName()

      if name == "rot_fast_ctrl":
        self._axis_rots[0] = obj.GetValue()
      elif name == "rot_slow_ctrl":
        self._axis_rots[1] = obj.GetValue()

      self.DrawRings()

    def OnSpinCell(self, event):
        obj = event.EventObject
        name = obj.GetName()

        self._cell[self._cell_control_names.index(name)] = obj.GetValue()

        self.DrawRings()

    def OnSpin(self, event):
        self.DrawRings()

    def OnSpaceGroup(self, event):
        obj = event.EventObject
        self._spacegroup = obj.GetValue()

        self.DrawRings()

    def OnClear(self, event):
        self.__del__()

    def _draw_rings_layer(self, dc, data, map_rel):
        """Draw a points layer.

        dc       the device context to draw on
        data     an iterable of point tuples:
                 (x, y, place, radius, colour, x_off, y_off, pdata)
        map_rel  points relative to map if True, MUST BE TRUE for lightweight
        Assumes all points are the same colour, saving 100's of ms.
        """

        assert map_rel is True
        if len(data) == 0:
            return
        (lon, lat, place, radius, colour, x_off, y_off, pdata) = data[0]

        scale = 2**self._pyslip.tiles.zoom_level

        # Draw points on map/view, using transparency if implemented.
        try:
            dc = wx.GCDC(dc)
        except NotImplementedError:
            pass
        dc.SetPen(wx.Pen(colour))
        dc.SetBrush(wx.Brush(colour, wx.TRANSPARENT))
        for (lon, lat, place, radius, colour, x_off, y_off, pdata) in data:
            (x, y) = self._pyslip.ConvertGeo2View((lon, lat))
            dc.DrawCircle(x, y, radius * scale)

    def DrawRings(self):
        frame = self.GetParent().GetParent()

        try:
            uc = symmetry(
                unit_cell=self._cell, space_group_symbol=str(self._spacegroup)
            )
            hkl_list = cctbx.miller.build_set(
                uc, False, d_min=self.d_min_ctrl.GetValue()
            )
        except Exception as e:
            frame.update_statusbar(str(e))
            return

        if self._show_hkl:
            hkl_list = hkl_list.common_set(
                cctbx.miller.set(
                    crystal_symmetry=uc, indices=self._show_hkl, anomalous_flag=False
                )
            )

        frame.update_statusbar(
            "%d %d %d %d %d %d, " % tuple(self._cell)
            + f"number of indices: {len(hkl_list.indices())}"
        )

        spacings = sorted(hkl_list.d_spacings(), key=lambda s: s[1], reverse=True)
        print(f"Printing spacings, len: {len(spacings)}")

        for d in spacings:
            print(d)
        d_values = [d[1] for d in spacings]

        detector = self._pyslip.tiles.raw_image.get_detector()
        #detector = copy.deepcopy(self._detector)

        pixel_size = detector[0].get_pixel_size()[
            0
        ]  # FIXME assumes square pixels, and that all panels use same pixel size
        origin_shift = (
            -1 * self._center[0] * pixel_size,
            self._center[1] * pixel_size,
            -1 * self.distance_ctrl.GetValue()
            )

        rot_fast_mx = col(self.fast_cache).axis_and_angle_as_r3_rotation_matrix(
            self._axis_rots[0], deg=True
            )
        rot_slow_mx = col(self.slow_cache).axis_and_angle_as_r3_rotation_matrix(
            self._axis_rots[1], deg=True
            )
        rot_mx = rot_fast_mx * rot_slow_mx
        new_origin = tuple(
            (a+b for a,b in zip(self.origin_cache, origin_shift))
            )

        detector.hierarchy().set_frame(
            rot_mx * col(self.fast_cache),
            rot_mx * col(self.slow_cache),
            new_origin
            )
        with open('modified.expt', 'w') as f:
            f.write(json.dumps({"detector": detector.to_dict()}, indent=2))

        wavelength = float(self.wavelength_ctrl.GetValue())
        distance = float(self.distance_ctrl.GetValue()) #MOVE DET

        frame.draw_resolution_rings(d_values=d_values)
        frame.pyslip.Update()
        panel = detector[0]
        fast = col(panel.get_fast_axis())
        slow = col(panel.get_slow_axis())
        norm = col(panel.get_normal())
        x = -panel.pixel_to_millimeter(self._center)[0]
        y = -panel.pixel_to_millimeter(self._center)[1]
        z = -(panel.get_distance() - distance)
        origin = (fast * x + slow * y + norm * z) + col(panel.get_origin())
        self.origin.SetLabel("Panel 0 origin: %f, %f, %f" % origin.elems)

