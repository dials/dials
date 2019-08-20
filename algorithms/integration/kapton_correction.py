from __future__ import absolute_import, division, print_function

import math

from dials.algorithms.shoebox import MaskCode
from dials.array_family import flex
from iotbx import phil
from scitbx import matrix

absorption_defs = """
  absorption_correction
    .multiple = True {
    apply = False
      .type = bool
      .help = must be supplied as a user-defined function with a specific interface (not documented)
    algorithm = fuller_kapton kapton_2019 other
      .type = choice
      .help = a specific absorption correction, or implementation thereof \
              kapton_2019 is a more general implementation of fuller_kapton \
              for use on single/multi-panel detectors
    fuller_kapton {
      xtal_height_above_kapton_mm {
        value = 0.02
          .type = float
          .help = height of the beam (or the irradiated crystal) above the kapton tape
        sigma = 0.01
          .type = float
        }
      rotation_angle_deg {
        value = 1.15
          .type = float
          .help = angle of the tape from vertical
        sigma = 0.1
          .type = float
        }
      kapton_half_width_mm {
        value = 1.5875
          .type = float
          .help = forward distance from irradiated crystal to edge of tape nearest detector
        sigma = 0.5
          .type = float
        }
      kapton_thickness_mm {
        value = 0.05
          .type = float
          .help = tape thickness
        sigma = 0.005
          .type = float
        }
      smart_sigmas = False
        .type = bool
        .help = apply spot-specific sigma corrections using kapton param sigmas
      within_spot_sigmas = True
        .type = bool
        .help = calculate initial per-spot sigmas based on variance across pixels in the spot.
        .help = turn this off to get a major speed-up
    }
  }"""

absorption_phil_scope = phil.parse(absorption_defs, process_includes=True)


class get_absorption_correction(object):
    def __init__(self):
        # Kapton, or polyimide.  C22H10N2O5 Density=1.43, Angle=90.deg
        # Photon Energy (eV), Atten Length (microns)
        data = """6000.00  482.643
   6070.00  500.286
   6140.00  518.362
   6210.00  536.896
   6280.00  555.873
   6350.00  575.302
   6420.00  595.191
   6490.00  615.552
   6560.00  636.382
   6630.00  657.691
   6700.00  679.484
   6770.00  701.758
   6840.00  724.521
   6910.00  747.791
   6980.00  771.561
   7050.00  795.846
   7120.00  820.646
   7190.00  845.963
   7260.00  871.812
   7330.00  898.183
   7400.00  925.082
   7470.00  952.535
   7540.00  980.535
   7610.00  1009.09
   7680.00  1038.18
   7750.00  1067.85
   7820.00  1098.08
   7890.00  1128.88
   7960.00  1160.25
   8030.00  1192.20
   8100.00  1224.76
   8170.00  1257.91
   8240.00  1291.67
   8310.00  1326.01
   8380.00  1360.98
   8450.00  1396.54
   8520.00  1432.72
   8590.00  1469.51
   8660.00  1506.93
   8730.00  1544.96
   8800.00  1583.65
   8870.00  1622.95
   8940.00  1662.90
   9010.00  1703.49
   9080.00  1744.72
   9150.00  1786.59
   9220.00  1829.13
   9290.00  1872.31
   9360.00  1916.16
   9430.00  1960.65
   9500.00  2005.82
   9570.00  2051.65
   9640.00  2098.16
   9710.00  2145.36
   9780.00  2193.22
   9850.00  2241.75
   9920.00  2290.95
   9990.00  2340.86
   10060.0  2391.49
   10130.0  2442.84
   10200.0  2494.86
   10270.0  2547.59
   10340.0  2601.02
   10410.0  2655.14
   10480.0  2709.98
   10550.0  2765.49
   10620.0  2821.73
   10690.0  2878.68
   10760.0  2936.31
   10830.0  2994.67
   10900.0  3053.72
   10970.0  3113.49
   11040.0  3173.96
   11110.0  3235.14
   11180.0  3297.07
   11250.0  3359.67
   11320.0  3423.01
   11390.0  3487.04
   11460.0  3551.76
   11530.0  3617.23
   11600.0  3683.38
   11670.0  3750.23
   11740.0  3817.81
   11810.0  3886.07
   11880.0  3955.05
   11950.0  4024.75
   12020.0  4095.11
   12090.0  4166.20
   12160.0  4237.96
   12230.0  4310.40
   12300.0  4383.60
   12370.0  4457.48
   12440.0  4532.02
   12510.0  4607.25
   12580.0  4683.14
   12650.0  4759.73
   12720.0  4837.01
   12790.0  4914.94
   12860.0  4993.54
   12930.0  5072.79
   13000.0  5152.69"""
        self.energy = flex.double()
        self.microns = flex.double()
        for line in data.split("\n"):
            tokens = line.strip().split()
            self.energy.append(float(tokens[0]))
            self.microns.append(float(tokens[1]))

    def __call__(self, wavelength_ang):
        # calculate energy in eV 12398.425 eV/Ang
        energy_eV = 12398.425 / wavelength_ang
        # interpolate the Henke tables downloaded from lbl.gov
        index_float = (
            (len(self.energy) - 1)
            * (energy_eV - self.energy[0])
            / (self.energy[-1] - self.energy[0])
        )
        fraction, int_idx = math.modf(index_float)
        int_idx = int(int_idx)

        microns = self.microns[int_idx] + fraction * (
            self.microns[int_idx + 1] - self.microns[int_idx]
        )
        return microns / 1000.0


class KaptonAbsorption(object):
    def __init__(
        self,
        height_mm,
        thickness_mm,
        half_width_mm,
        rotation_angle_deg,
        detector_dist_mm=None,
        pixel_size_mm=None,
        wavelength_ang=None,
        size_fast=None,
        size_slow=None,
        detector=None,
        beam=None,
    ):
        self.height_mm = height_mm  # tool controlled
        self.thickness_mm = thickness_mm  # tool controlled
        self.half_width_mm = half_width_mm  # tool controlled
        self.angle_rad = rotation_angle_deg * math.pi / 180.0  # tool controlled
        self.detector_dist_mm = detector_dist_mm
        self.pixel_size_mm = pixel_size_mm
        self.wavelength_ang = wavelength_ang
        self.size_fast = size_fast
        self.size_slow = size_slow
        self.detector = detector
        self.beam = beam

        if (
            False and self.detector is not None
        ):  ### XXX not implemented -- impose fast and slow axes for now
            self.fast = matrix.col(self.detector[0].get_fast_axis())
            self.slow = matrix.col(self.detector[0].get_slow_axis())
            print("fast and slow detector axes:", tuple(self.fast), tuple(self.slow))
        else:
            self.fast = matrix.col((1, 0, 0))
            self.slow = matrix.col((0, -1, 0))
        if self.beam is not None:
            self.beam_direction = matrix.col(self.beam.get_s0()).normalize()
        else:
            self.beam_direction = matrix.col((0, 0, -1))

        # determine absorption coeff (mm-1) through kapton for a given X-ray energy
        G = get_absorption_correction()
        attenuation_length_mm = G(self.wavelength_ang)
        self.abs_coeff = 1 / attenuation_length_mm

        # determine zones of different absorption behavior
        self.surface_normal = matrix.col((-1, 0, 0)).rotate_around_origin(
            axis=self.beam_direction, angle=-self.angle_rad
        )

        # ray tracing
        self.edge_of_tape_normal = self.beam_direction
        self.surface1_point_mm = -self.surface_normal * self.height_mm
        self.surface2_point_mm = -self.surface_normal * (
            self.height_mm + self.thickness_mm
        )
        self.surface3_point_mm = (
            self.surface1_point_mm + self.edge_of_tape_normal * self.half_width_mm
        )
        self.sn1 = self.surface1_point_mm.dot(self.surface_normal)
        self.sn2 = self.surface2_point_mm.dot(self.surface_normal)
        self.sn3 = self.surface3_point_mm.dot(self.edge_of_tape_normal)

    def abs_correction(self, s1):
        try:
            # let's get the unit vector along the direction of the X-ray
            dsurf1 = self.sn1 / (s1.dot(self.surface_normal))
            dsurf2 = self.sn2 / (s1.dot(self.surface_normal))
            dsurf3 = self.sn3 / (s1.dot(self.edge_of_tape_normal))

            # determine path length through kapton tape
            if dsurf3 < dsurf1 or dsurf1 < 0:
                kapton_path_mm = 0.0
            elif dsurf3 < dsurf2:
                kapton_path_mm = dsurf3 - dsurf1
            else:
                kapton_path_mm = dsurf2 - dsurf1

            # determine absorption correction
            absorption_correction = 1 / math.exp(
                -self.abs_coeff * kapton_path_mm
            )  # unitless, >=1
            return absorption_correction
        except ZeroDivisionError:
            return 0

    def abs_correction_flex(self, s1_flex):
        s1_flex_length = len(s1_flex)
        surface_normal = flex.vec3_double(s1_flex_length, self.surface_normal)
        edge_of_tape_normal = flex.vec3_double(s1_flex_length, self.edge_of_tape_normal)

        dsurf1 = flex.double(len(s1_flex), 0)
        dsurf2 = flex.double(len(s1_flex), 0)
        dot_product = s1_flex.dot(surface_normal)
        sel = dot_product != 0
        dsurf1.set_selected(sel, self.sn1 / dot_product.select(sel))
        dsurf2.set_selected(sel, self.sn2 / dot_product.select(sel))
        dsurf3 = self.sn3 / (s1_flex.dot(edge_of_tape_normal))

        # determine path length through kapton tape
        kapton_path_mm = flex.double(s1_flex_length, 0)
        unshadowed_sel = (dsurf3 < dsurf1) | (dsurf1 < 0)
        nearsel = ~unshadowed_sel & (dsurf3 < dsurf2)
        kapton_path_mm.set_selected(nearsel, (dsurf3 - dsurf1).select(nearsel))
        farsel = ~unshadowed_sel & (dsurf3 >= dsurf2)
        kapton_path_mm.set_selected(farsel, (dsurf2 - dsurf1).select(farsel))

        # determine absorption correction
        absorption_correction = 1 / flex.exp(
            -self.abs_coeff * kapton_path_mm
        )  # unitless, >=1
        return absorption_correction

    def abs_bounding_lines(self, as_xy_ints=False):
        if not hasattr(
            self, "segments"
        ):  # avoid recomputing the bounding lines on subsequent calls
            # determine the beam center in mm
            self.beam_ctr_mm_on_detector = matrix.col(
                (
                    self.detector[0].get_beam_centre_px(self.beam_direction)[1]
                    * self.pixel_size_mm,
                    self.detector[0].get_beam_centre_px(self.beam_direction)[0]
                    * self.pixel_size_mm,
                    0,
                )
            )
            # dials format coordinates in pixels of the detector corner
            self.det_zero_mm = (
                self.detector_dist_mm * self.beam_direction
                - self.beam_ctr_mm_on_detector
            )
            # calculate the positions of the boundaries between zones with different absorption behavior
            # get a vector orthogonal to the surface normal and edge normal (along the direction the tape travels)
            self.tape_dir = self.surface_normal.cross(self.edge_of_tape_normal)
            # and one normal to the detector surface
            self.det_normal = self.slow.cross(self.fast)
            # get one more point on the surface of the tape
            self.surface4_point_mm = (
                self.surface3_point_mm - self.surface_normal * self.thickness_mm
            )
            # get the unit vectors normal to the planes of the absorption minimum and maximum
            self.min_normal = self.tape_dir.cross(
                self.surface3_point_mm
            )  # i.e. normal to plane through near edge of tape
            self.max_normal = self.tape_dir.cross(
                self.surface4_point_mm
            )  # i.e. normal to plane through far edge of tape
            self.center_normal_1 = self.tape_dir.cross(
                self.edge_of_tape_normal
            )  # crosshairs part 1
            self.center_normal_2 = self.surface_normal.cross(
                self.edge_of_tape_normal
            )  # crosshairs part 2
            # get some intermediate values (see LaTeX documentation)
            min_det_dot = self.min_normal - self.det_normal
            max_det_dot = self.max_normal - self.det_normal
            # ctr1_det_dot = self.center_normal_1 - self.det_normal # crosshairs part 1
            # ctr2_det_dot = self.center_normal_2 - self.det_normal # crosshairs part 2
            min_max_det_add = self.det_zero_mm.dot(self.det_normal)
            # get the distances along the detector edges to the minimum and maximum absorption edges
            def get_intersection(offset, direction, min_or_max_det_dot, dimension):
                # offset is a number of pixels from the detector origin
                # direction is the fast or slow axis of the detector
                # min_or_max_det_dot is a quantity to be used in a dot product (see documentation)
                # dimension is the number of pixels in the panel along the direction specified
                try:
                    num = (
                        -(self.det_zero_mm + offset * self.pixel_size_mm).dot(
                            min_or_max_det_dot
                        )
                        - min_max_det_add
                    )
                    denom = direction.dot(min_or_max_det_dot)
                    dist = num / denom
                    if (
                        0 <= dist and dist <= dimension * self.pixel_size_mm
                    ):  # point is on a detector edge
                        return (
                            self.det_zero_mm
                            + offset * self.pixel_size_mm
                            + dist * direction
                        )
                    else:
                        return None
                except ZeroDivisionError:
                    return None

            zero = matrix.col((0, 0, 0))
            self.segments = []
            for edge in [min_det_dot, max_det_dot]:
                # for edge in [min_det_dot, max_det_dot, ctr1_det_dot, ctr2_det_dot]: # crosshairs
                # calculate all the possible intersections on the detector edges
                d_top_fast = get_intersection(zero, self.fast, edge, self.size_fast)
                d_bottom_fast = get_intersection(
                    self.size_slow * self.slow, self.fast, edge, self.size_fast
                )
                d_left_slow = get_intersection(zero, self.slow, edge, self.size_slow)
                d_right_slow = get_intersection(
                    self.size_fast * self.fast, self.slow, edge, self.size_slow
                )
                # collect the ones that aren't None
                points = set()
                for p in [d_top_fast, d_bottom_fast, d_left_slow, d_right_slow]:
                    if p is not None:
                        points.add(p)
                assert len(points) <= 2, (
                    "There seems to be a plane intersection at a corner of the detector "
                    + "that is registering as an intersection on both the fast and the slow axes with slightly different "
                    + "coordinates. Please consider figuring out how to discard one of the two."
                )
                assert len(points) == 2, "The edges of the detector were not found."
                self.segments.append(tuple(points))
        if as_xy_ints:
            if not hasattr(self, "segments_as_xy"):
                self.segments_as_xy = []
                for seg in self.segments:
                    seg_as_xy = []
                    for point in seg:
                        as_xy = point - self.det_zero_mm
                        x_int = int(as_xy.dot(self.slow) / self.pixel_size_mm)
                        y_int = int(as_xy.dot(self.fast) / self.pixel_size_mm)
                        seg_as_xy += [x_int, y_int]
                    self.segments_as_xy.append(tuple(seg_as_xy))
            return self.segments_as_xy
        else:
            return self.segments


class image_kapton_correction(object):
    def __init__(
        self,
        panel_size_px=None,
        pixel_size_mm=None,
        detector_dist_mm=None,
        wavelength_ang=None,
        reflections_sele=None,
        params=None,
        expt=None,
        refl=None,
        smart_sigmas=True,
        logger=None,
    ):
        self.panel_size_px = panel_size_px
        self.pixel_size_mm = pixel_size_mm
        self.detector_dist_mm = detector_dist_mm
        self.wavelength_ang = wavelength_ang
        self.reflections_sele = reflections_sele
        self.params = params
        self.expt = expt
        self.refl = refl
        self.smart_sigmas = smart_sigmas
        self.logger = logger
        self.extract_params()

    def extract_params(self):
        h = self.params.xtal_height_above_kapton_mm.value
        t = self.params.kapton_thickness_mm.value
        w = self.params.kapton_half_width_mm.value
        a = self.params.rotation_angle_deg.value
        self.kapton_params = (h, t, w, a)

        if self.smart_sigmas:
            sig_h = self.params.xtal_height_above_kapton_mm.sigma
            sig_t = self.params.kapton_thickness_mm.sigma
            sig_w = self.params.kapton_half_width_mm.sigma
            sig_a = self.params.rotation_angle_deg.sigma
            self.kapton_params_sigmas = (sig_h, sig_t, sig_w, sig_a)
            assert all(
                sig >= 0 for sig in self.kapton_params_sigmas
            ), "Kapton param sigmas must be nonnegative"
            self.kapton_params_maxes = [
                [
                    self.kapton_params[i] + self.kapton_params_sigmas[j]
                    if j == i
                    else self.kapton_params[i]
                    for i in range(4)
                ]
                for j in range(4)
            ]
            self.kapton_params_mins = [
                [
                    max(self.kapton_params[i] - self.kapton_params_sigmas[j], 0.001)
                    if j == i
                    else self.kapton_params[i]
                    for i in range(3)
                ]
                + [a]
                for j in range(3)
            ] + [[self.kapton_params[i] for i in range(3)] + [a - sig_a]]

    def __call__(self, plot=False):
        def correction_and_within_spot_sigma(params_version, variance_within_spot=True):
            absorption = KaptonAbsorption(
                params_version[0],
                params_version[1],
                params_version[2],
                params_version[3],
                self.detector_dist_mm,
                self.pixel_size_mm,
                self.wavelength_ang,
                *map(float, self.panel_size_px)
            )
            detector = self.expt.detector
            beam = self.expt.beam
            # y_max = int(detector[0].millimeter_to_pixel(detector[0].get_image_size())[1])
            s0_fast, s0_slow = map(int, detector[0].get_beam_centre_px(beam.get_s0()))

            absorption_corrections = flex.double()
            absorption_sigmas = (
                flex.double()
            )  # std dev of corrections for pixels within a spot, default sigma

            if variance_within_spot:
                mask_code = MaskCode.Foreground | MaskCode.Valid
                for iref in range(len(self.reflections_sele)):
                    kapton_correction_vector = flex.double()
                    # foreground: integration mask
                    shoebox = self.reflections_sele[iref]["shoebox"]
                    foreground = (
                        (shoebox.mask.as_1d() & mask_code) == mask_code
                    ).iselection()
                    width = shoebox.xsize()
                    fast_coords = (foreground % width).as_int()  # within spot
                    slow_coords = (foreground / width).as_int()  # within spot
                    f_absolute = fast_coords + shoebox.bbox[0]  # relative to detector
                    s_absolute = slow_coords + shoebox.bbox[2]  # relative to detector
                    lab_coords = detector[0].get_lab_coord(
                        detector[0].pixel_to_millimeter(
                            flex.vec2_double(
                                f_absolute.as_double(), s_absolute.as_double()
                            )
                        )
                    )
                    s1 = lab_coords.each_normalize()
                    kapton_correction_vector.extend(absorption.abs_correction_flex(s1))
                    average_kapton_correction = flex.mean(kapton_correction_vector)
                    absorption_corrections.append(average_kapton_correction)
                    try:
                        spot_px_stddev = flex.mean_and_variance(
                            kapton_correction_vector
                        ).unweighted_sample_standard_deviation()
                    except Exception:
                        assert (
                            len(kapton_correction_vector) == 1
                        ), "stddev could not be calculated"
                        spot_px_stddev = 0
                    absorption_sigmas.append(spot_px_stddev)
                return absorption_corrections, absorption_sigmas
            else:
                s1_flex = self.reflections_sele["s1"].each_normalize()
                absorption_corrections = absorption.abs_correction_flex(s1_flex)
                return absorption_corrections, None

        # loop through modified Kapton parameters to get alternative corrections and estimate sigmas as
        # maximum variation between these versions of the corrections, on a per-spot basis, or the standard
        # deviation within a single spot, whichever is larger.
        self.logger.info("Calculating kapton corrections to integrated intensities...")
        corrections, sigmas = correction_and_within_spot_sigma(
            self.kapton_params, variance_within_spot=self.params.within_spot_sigmas
        )
        if self.smart_sigmas:
            for p in self.kapton_params_mins + self.kapton_params_maxes:
                self.logger.info("Calculating smart sigmas...")
                modif_corrections, _ = correction_and_within_spot_sigma(
                    p, variance_within_spot=False
                )
                perturbed = flex.abs(corrections - modif_corrections)
                if sigmas is None:
                    sigmas = perturbed
                else:
                    replace_sel = perturbed > sigmas
                    sigmas.set_selected(replace_sel, perturbed.select(replace_sel))
        if plot:
            from matplotlib import pyplot as plt

            for (title, data) in [("corrections", corrections), ("sigmas", sigmas)]:
                n, bins, patches = plt.hist(data, 20)
                plt.title(title)
                plt.show()
        if self.logger is not None:
            self.logger.info(
                "Returning absorption corrections and sigmas for %d spots"
                % len(corrections)
            )
        return corrections, sigmas


class multi_kapton_correction(object):
    def __init__(self, experiments, integrated, kapton_params, logger=None):
        self.experiments = experiments
        self.reflections = integrated
        self.params = kapton_params
        self.logger = logger

    def __call__(self):
        self.corrected_reflections = flex.reflection_table()
        for expt, refl in zip(
            self.experiments, self.reflections.split_by_experiment_id()
        ):
            # extract experiment details
            detector = expt.detector
            assert len(detector) == 1  # for now
            panel = detector[0]
            panel_size_px = panel.get_image_size()
            pixel_size_mm = panel.get_pixel_size()[0]
            detector_dist_mm = panel.get_distance()
            beam = expt.beam
            wavelength_ang = beam.get_wavelength()

            # exclude reflections with no foreground pixels
            refl_valid = refl.select(
                refl["num_pixels.valid"] > 0 and refl["num_pixels.foreground"] > 0
            )
            refl_zero = refl_valid.select(refl_valid["intensity.sum.value"] == 0)
            refl_nonzero = refl_valid.select(refl_valid["intensity.sum.value"] != 0)

            def correct(refl_sele, smart_sigmas=True):
                kapton_correction = image_kapton_correction(
                    panel_size_px=panel_size_px,
                    pixel_size_mm=pixel_size_mm,
                    detector_dist_mm=detector_dist_mm,
                    wavelength_ang=wavelength_ang,
                    reflections_sele=refl_sele,
                    params=self.params,
                    expt=expt,
                    refl=refl,
                    smart_sigmas=smart_sigmas,
                    logger=self.logger,
                )

                k_corr, k_sigmas = kapton_correction()
                refl_sele["kapton_absorption_correction"] = k_corr
                if smart_sigmas:
                    refl_sele["kapton_absorption_correction_sigmas"] = k_sigmas
                    # apply corrections and propagate error
                    # term1 = (sig(C)/C)^2
                    # term2 = (sig(Imeas)/Imeas)^2
                    # I' = C*I
                    # sig^2(I') = (I')^2*(term1 + term2)
                    integrated_data = refl_sele["intensity.sum.value"]
                    integrated_variance = refl_sele["intensity.sum.variance"]
                    integrated_sigma = flex.sqrt(integrated_variance)
                    term1 = flex.pow(k_sigmas / k_corr, 2)
                    term2 = flex.pow(integrated_sigma / integrated_data, 2)
                    integrated_data *= k_corr
                    integrated_variance = flex.pow(integrated_data, 2) * (term1 + term2)
                    refl_sele["intensity.sum.value"] = integrated_data
                    refl_sele["intensity.sum.variance"] = integrated_variance
                    # order is purposeful: the two lines above require that integrated_data
                    # has already been corrected!
                else:
                    refl_sele["intensity.sum.value"] *= k_corr
                    refl_sele["intensity.sum.variance"] *= (k_corr) ** 2
                return refl_sele

            if len(refl_zero) > 0 and self.params.smart_sigmas:
                # process nonzero intensity reflections with smart sigmas as requested
                # but turn them off for zero intensity reflections to avoid a division by zero
                # during error propogation. Not at all certain this is the best way.
                self.corrected_reflections.extend(
                    correct(refl_nonzero, smart_sigmas=True)
                )
                self.corrected_reflections.extend(
                    correct(refl_zero, smart_sigmas=False)
                )
            else:
                self.corrected_reflections.extend(
                    correct(refl_valid, smart_sigmas=self.params.smart_sigmas)
                )

        return self.experiments, self.corrected_reflections
