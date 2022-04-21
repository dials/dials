"""Simulate a rotation dataset with a smoothly-varying beam position for
refinement testing. Script based on tst_nanoBragg_basic.py"""


from __future__ import annotations

import math

from cctbx.eltbx import henke
from dxtbx.model import (
    BeamFactory,
    Crystal,
    DetectorFactory,
    GoniometerFactory,
    ScanFactory,
)
from iotbx import pdb
from scitbx import matrix

_pdb_lines = """HEADER TEST
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""


class Simulation:
    def __init__(self, override_fdp=None):

        # Set up detector
        distance = 100
        pixel_size = 0.1
        image_size = (1000, 1000)
        beam_centre_mm = (
            pixel_size * image_size[0] / 2,
            pixel_size * image_size[1] / 2,
        )
        self.detector = DetectorFactory().simple(
            "CCD",
            distance,
            beam_centre_mm,
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
        )

        # Set up beam
        self.beam = BeamFactory().simple(wavelength=1)

        # Set up scan
        sequence_width = 90.0
        osc_start = 0.0
        image_width = 0.2
        oscillation = (osc_start, image_width)

        nframes = int(math.ceil(sequence_width / image_width))
        image_range = (1, nframes)
        exposure_times = 0.0
        epochs = [0] * nframes
        self.scan = ScanFactory().make_scan(
            image_range, exposure_times, oscillation, epochs, deg=True
        )

        # Set up goniometer
        self.goniometer = GoniometerFactory.known_axis(self.detector[0].get_fast_axis())

        # Set up simulated structure factors
        self.sfall = self.fcalc_from_pdb(
            resolution=1.6, algorithm="direct", override_fdp=override_fdp
        )

        # Set up crystal
        self.crystal = Crystal(
            real_space_a=(50, 0, 0),
            real_space_b=(0, 60, 0),
            real_space_c=(0, 0, 70),
            space_group_symbol="P1",
        )
        axis = matrix.col(
            elems=(-0.14480368275412925, -0.6202131724405818, -0.7709523423610766)
        )
        self.crystal.set_U(
            axis.axis_and_angle_as_r3_rotation_matrix(angle=0.625126343998969)
        )

    def fcalc_from_pdb(self, resolution, algorithm=None, override_fdp=None):
        pdb_inp = pdb.input(source_info=None, lines=_pdb_lines)
        xray_structure = pdb_inp.xray_structure_simple()
        wavelength = self.beam.get_wavelength()
        #
        # take a detour to calculate anomalous contribution of every atom
        scatterers = xray_structure.scatterers()
        for sc in scatterers:
            expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
            sc.fp = expected_henke.fp()
            sc.fdp = override_fdp if override_fdp is not None else expected_henke.fdp()

        # how do we do bulk solvent?
        primitive_xray_structure = xray_structure.primitive_setting()
        P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
        fcalc = P1_primitive_xray_structure.structure_factors(
            d_min=resolution, anomalous_flag=True, algorithm=algorithm
        ).f_calc()
        return fcalc.amplitudes()

    def set_varying_beam(self, along="fast", npixels_drift=5):
        assert along in ["fast", "slow", "both"]
        num_scan_points = self.scan.get_num_images() + 1
        s0 = matrix.col(self.beam.get_s0())
        beam_centre_px = self.detector[0].get_beam_centre_px(s0)
        if along == "fast":
            start_beam_centre = (
                beam_centre_px[0] - npixels_drift / 2,
                beam_centre_px[1],
            )
            end_beam_centre = (beam_centre_px[0] + npixels_drift / 2, beam_centre_px[1])
        elif along == "slow":
            start_beam_centre = (
                beam_centre_px[0],
                beam_centre_px[1] - npixels_drift / 2,
            )
            end_beam_centre = (beam_centre_px[0], beam_centre_px[1] + npixels_drift / 2)
        elif along == "both":
            offset = math.sqrt(2.0) * npixels_drift / 4.0
            start_beam_centre = (beam_centre_px[0] - offset, beam_centre_px[1] - offset)
            end_beam_centre = (beam_centre_px[0] + offset, beam_centre_px[1] + offset)

        start_lab = matrix.col(self.detector[0].get_pixel_lab_coord(start_beam_centre))
        end_lab = matrix.col(self.detector[0].get_pixel_lab_coord(end_beam_centre))
        axis = start_lab.cross(end_lab).normalize()
        full_angle = start_lab.angle(end_lab)
        angle_step = full_angle / self.scan.get_num_images()
        angles = [e * angle_step for e in range(num_scan_points)]

        start_s0 = start_lab.normalize() * s0.length()
        s0_list = [start_s0.rotate_around_origin(axis=axis, angle=e) for e in angles]

        self.beam.set_s0_at_scan_points(s0_list)
