from __future__ import division, print_function, absolute_import

import logging
import math

from dials.array_family import flex
from dials.algorithms.integration.kapton_correction import get_absorption_correction
from dials.algorithms.shoebox import MaskCode
from scitbx.matrix import col

logging.basicConfig()
logger = logging.getLogger(__name__)


class KaptonTape_2019(object):
    """Class for defining Kapton tape using dxtbx models and finding the path through the tape traversed by s1 vector"""

    def __init__(
        self,
        height_mm,
        thickness_mm,
        half_width_mm,
        rotation_angle_deg,
        wavelength_ang=None,
    ):
        self.height_mm = h = height_mm  # plugin controlled
        self.thickness_mm = t = thickness_mm  # plugin controlled
        self.half_width_mm = w = half_width_mm  # plugin controlled
        self.tape_depth_mm = a = half_width_mm * 20.0  # Big number in mm
        self.angle_rad = rotation_angle_deg * math.pi / 180.0  # plugin controlled
        self.wavelength_ang = wavelength_ang
        self.num_pixels = 5000  # number of pixels to put in fictitious kapton faces
        # Now set up the kapton physical model using dxtbx detector objects
        #
        # determine absorption coeff (mm-1) through kapton for a given X-ray energy
        G = get_absorption_correction()
        attenuation_length_mm = G(self.wavelength_ang)
        self.abs_coeff = 1 / attenuation_length_mm

        def create_kapton_face(ori, fast, slow, image_size, pixel_size, name):
            """ Create a face of the kapton as a dxtbx detector object """
            from dxtbx.model import Detector

            d = Detector()
            p = d.add_panel()
            p.set_local_frame(fast.elems, slow.elems, ori.elems)
            p.set_pixel_size((pixel_size, pixel_size))
            p.set_image_size(image_size)
            p.set_trusted_range((-1, 2e6))
            p.set_name("KAPTON_%s" % name)
            return d

        # Set up the bounding box of the kapton
        rA = col((t / 2, a / 2, -w))
        rB = col((t / 2, a / 2, w))
        rC = col((t / 2, -a / 2, w))
        rD = col((t / 2, -a / 2, -w))
        rE = col((-t / 2, a / 2, -w))
        rF = col((-t / 2, a / 2, w))
        rG = col((-t / 2, -a / 2, w))
        rH = col((-t / 2, -a / 2, -w))
        # Now add a rotation
        # rotation about z-axis
        rot_axis = col((0, 0, 1))
        rot_mat = rot_axis.axis_and_angle_as_r3_rotation_matrix(
            self.angle_rad, deg=False
        )

        rA = rot_mat * rA
        rB = rot_mat * rB
        rC = rot_mat * rC
        rD = rot_mat * rD
        rE = rot_mat * rE
        rF = rot_mat * rF
        rG = rot_mat * rG
        rH = rot_mat * rH
        # Now add an offset to all the points in the direction normal to the ABCD plane
        rAB = rA - rB
        rCB = rC - rB
        # Get normal in the direction of the drop
        drop_normal = rCB.cross(rAB).normalize()
        offset = -(h + t / 2) * drop_normal
        rA = rA + offset
        rB = rB + offset
        rC = rC + offset
        rD = rD + offset
        rE = rE + offset
        rF = rF + offset
        rG = rG + offset
        rH = rH + offset
        # Store these edge points if needed later
        self.edge_points = [rA, rB, rC, rD, rE, rF, rG, rH]
        # Now set up the 6 faces
        faces = []
        px_off = 0
        fast = (rA - rB).normalize()
        slow = (rC - rB).normalize()
        pixel_size = self.half_width_mm * 2.0 / self.num_pixels
        image_size = (
            self.num_pixels + px_off,
            int(self.tape_depth_mm / pixel_size) + px_off,
        )
        # Face YZ0
        ori = rB
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "yz0"))
        # face YZ1
        ori = rF
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "yz1"))
        #
        fast = (rB - rF).normalize()
        slow = (rG - rF).normalize()
        pixel_size = self.thickness_mm / self.num_pixels
        image_size = (
            self.num_pixels + px_off,
            int(self.tape_depth_mm / pixel_size) + px_off,
        )
        # Face XY0
        ori = rF
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "xy0"))
        # face XY1
        ori = rE
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "xy1"))
        #
        fast = (rA - rE).normalize()
        slow = (rF - rE).normalize()
        pixel_size = self.thickness_mm / self.num_pixels
        image_size = (
            self.num_pixels + px_off,
            int(self.half_width_mm * 2.0 / pixel_size) + px_off,
        )
        # Face XZ0
        ori = rE
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "xz0"))
        # face XY1
        ori = rH
        faces.append(create_kapton_face(ori, fast, slow, image_size, pixel_size, "xz1"))
        #
        self.faces = faces

    def get_kapton_path_mm(self, s1):
        """ Get kapton path length traversed by an s1 vecto. If no kapton intersection or just touches the edge,
            then None is returned"""
        intersection_points = []
        # determine path length through kapton tape
        for face in self.faces:
            try:
                px = [int(p) for p in face[0].get_ray_intersection_px(s1)]
                # faces are single panel anyway
                if (
                    px[0] < 0
                    or px[1] < 0
                    or px[0] > face[0].get_image_size()[0]
                    or px[1] > face[0].get_image_size()[1]
                ):
                    continue
                intersection_points.append(
                    face[0].get_lab_coord(face[0].get_ray_intersection(s1))
                )  # faces are single panel anyway
            except RuntimeError:
                pass

        if not intersection_points:
            return 0.0
        if len(intersection_points) == 1:
            logger.warning(
                "Ray not intersecting 2 faces.Either s1 vector does not intersect kapton or just touches edge of kapton tape \
                    No correction is needed for those cases."
            )
            return 0.0

        n_intersection_points = len(intersection_points)
        kapton_path_mm = []
        for ii in range(n_intersection_points - 1):
            for jj in range(ii + 1, n_intersection_points):
                kapton_path_mm.append(
                    (
                        col(intersection_points[ii]) - col(intersection_points[jj])
                    ).length()
                )
        return max(kapton_path_mm)

    def abs_correction(self, s1):
        """ Compute absorption correction using beers law. Takes in a tuple for a single s1 vector and does absorption correction"""
        kapton_path_mm = self.get_kapton_path_mm(s1)
        if kapton_path_mm is not None:
            absorption_correction = 1 / math.exp(
                -self.abs_coeff * kapton_path_mm
            )  # unitless, >=1
        return absorption_correction

    def abs_correction_flex(self, s1_flex):
        """ Compute the absorption correction using beers law. Takes in a flex array of s1 vectors, determines path lengths for each
            and then determines absorption correction for each s1 vector """
        kapton_faces = self.faces
        from dials.algorithms.integration import get_kapton_path_cpp

        # new style, much faster
        kapton_path_mm = get_kapton_path_cpp(kapton_faces, s1_flex)
        # old style, really slow
        # for s1 in s1_flex:
        #  kapton_path_mm.append(self.get_kapton_path_mm(s1))
        # determine absorption correction
        if kapton_path_mm is not None:
            absorption_correction = 1 / flex.exp(
                -self.abs_coeff * kapton_path_mm
            )  # unitless, >=1
        return absorption_correction

    def distance_of_point_from_line(self, r0, r1, r2):
        """ Evaluates distance between point and a line between two points
            Note that implementation ignores z dimension"""
        x0, y0, z0 = r0
        x1, y1, z1 = r1
        x2, y2, z2 = r2
        num = (y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1
        denom = math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1))
        return abs(num) / denom

    def abs_bounding_lines_in_mm(self, detector):
        """ Return bounding lines of kapton """
        # first get bounding directions from detector:
        detz = flex.mean(flex.double([panel.get_origin()[2] for panel in detector]))
        edges = []
        for ii, panel in enumerate(detector):
            f_size, s_size = panel.get_image_size()
            for point in [(0, 0), (0, s_size), (f_size, 0), (f_size, s_size)]:
                x, y = panel.get_pixel_lab_coord(point)[0:2]
                edges.append((x, y, detz))
        # Use the idea that the corners of the detector are end points of the diagonal and will be the
        # top 2 max dimension among all end points
        dlist = flex.double()
        dlist_idx = []
        n_edges = len(edges)
        for ii in range(n_edges - 1):
            for jj in range(ii + 1, n_edges):
                pt_1 = col(edges[ii])
                pt_2 = col(edges[jj])
                distance = (pt_1 - pt_2).length()
                dlist.append(distance)
                dlist_idx.append((ii, jj))
        sorted_idx = flex.sort_permutation(dlist, reverse=True)

        edge_pts = [
            edges[dlist_idx[sorted_idx[0]][0]],
            edges[dlist_idx[sorted_idx[1]][0]],
            edges[dlist_idx[sorted_idx[0]][1]],
            edges[dlist_idx[sorted_idx[1]][1]],
        ]

        self.detector_edges = edge_pts
        # Now get the maximum extent of the intersection of the rays with the detector
        all_ints = []
        kapton_path_list = []

        for ii, edge_point in enumerate(self.edge_points):
            s1 = edge_point.normalize()
            kapton_path_mm = self.get_kapton_path_mm(s1)
            for panel in detector:
                try:
                    x_int, y_int = panel.get_lab_coord(panel.get_ray_intersection(s1))[
                        0:2
                    ]
                except RuntimeError:
                    pass
                int_point = (x_int, y_int, detz)
                # Arbitrary tolerance of couple of pixels otherwise these points were getting clustered together
                tolerance = min(panel.get_pixel_size()) * 2.0
                if (
                    sum(
                        [
                            (col(trial_pt) - col(int_point)).length() <= tolerance
                            for trial_pt in all_ints
                        ]
                    )
                    == 0
                ):
                    all_ints.append(int_point)
                kapton_path_list.append(kapton_path_mm)
        # Use the idea that the extreme edges of the intersection points are end points of the diagonal and will be the
        # top 2 max dimension among all end points
        dlist = flex.double()
        dlist_idx = []
        n_edges = len(all_ints)
        for ii in range(n_edges - 1):
            pt_1 = col(all_ints[ii])
            for jj in range(ii + 1, n_edges):
                pt_2 = col(all_ints[jj])
                distance = (pt_1 - pt_2).length()
                dlist.append(distance)
                dlist_idx.append((ii, jj))
        sorted_idx = flex.sort_permutation(dlist, reverse=True)

        int_edge_pts = [
            all_ints[dlist_idx[sorted_idx[0]][0]],
            all_ints[dlist_idx[sorted_idx[1]][0]],
            all_ints[dlist_idx[sorted_idx[0]][1]],
            all_ints[dlist_idx[sorted_idx[1]][1]],
        ]

        # Sort out the edge points and the int_edge_points which are on the same side
        kapton_edge_1 = (col(int_edge_pts[0]) - col(int_edge_pts[1])).normalize()
        kapton_edge_2 = (col(int_edge_pts[2]) - col(int_edge_pts[3])).normalize()
        min_loss_func = -999.9
        edge_idx = None
        for edge_idx_combo in [(0, 1, 2, 3), (0, 3, 1, 2)]:
            side_1 = (
                col(edge_pts[edge_idx_combo[0]]) - col(edge_pts[edge_idx_combo[1]])
            ).normalize()
            side_2 = (
                col(edge_pts[edge_idx_combo[2]]) - col(edge_pts[edge_idx_combo[3]])
            ).normalize()
            loss_func = abs(kapton_edge_1.dot(side_1)) + abs(kapton_edge_2.dot(side_2))
            if loss_func > min_loss_func:
                edge_idx = edge_idx_combo
                min_loss_func = loss_func
        # Make sure the edges of the detector and the kapton are in the same orientation
        # first for kapton edge 1
        side_1 = (col(edge_pts[edge_idx[0]]) - col(edge_pts[edge_idx[1]])).normalize()
        side_2 = (col(edge_pts[edge_idx[2]]) - col(edge_pts[edge_idx[3]])).normalize()
        v1 = kapton_edge_1.dot(side_1)
        v2 = kapton_edge_2.dot(side_2)
        if v1 < 0.0:
            edge_idx = (edge_idx[1], edge_idx[0], edge_idx[2], edge_idx[3])
        if v2 < 0.0:
            edge_idx = (edge_idx[0], edge_idx[1], edge_idx[3], edge_idx[2])

        # Now make sure the edges and the kapton lines are on the right side (i.e not swapped).
        # Let's look at edge_idx[0:2] i,e the first edge of detector parallel to the kapton
        pt1 = edge_pts[edge_idx[0]]
        pt2 = edge_pts[edge_idx[1]]
        # Now find the distance between each of these points and the kapton lines.
        d1_kapton_1 = self.distance_of_point_from_line(
            pt1, int_edge_pts[0], int_edge_pts[1]
        )
        d1_kapton_2 = self.distance_of_point_from_line(
            pt1, int_edge_pts[2], int_edge_pts[3]
        )
        d2_kapton_1 = self.distance_of_point_from_line(
            pt2, int_edge_pts[0], int_edge_pts[1]
        )
        d2_kapton_2 = self.distance_of_point_from_line(
            pt2, int_edge_pts[2], int_edge_pts[3]
        )
        if d1_kapton_1 < d1_kapton_2:  # closer to max than edge
            assert (
                d2_kapton_1 < d2_kapton_2
            ), "Distance mismatch. Edge of detector might be on wrong side of kapton tape ... please check"
            pair_values = [
                (
                    edge_pts[edge_idx[0]],
                    edge_pts[edge_idx[1]],
                    edge_pts[edge_idx[2]],
                    edge_pts[edge_idx[3]],
                ),
                (int_edge_pts[0], int_edge_pts[1], int_edge_pts[2], int_edge_pts[3]),
            ]
        else:
            pair_values = [
                (
                    edge_pts[edge_idx[0]],
                    edge_pts[edge_idx[1]],
                    edge_pts[edge_idx[2]],
                    edge_pts[edge_idx[3]],
                ),
                (int_edge_pts[2], int_edge_pts[3], int_edge_pts[0], int_edge_pts[1]),
            ]

        return pair_values

    def abs_bounding_lines_on_image(self, detector):
        pair_values = self.abs_bounding_lines_in_mm(detector)
        r0, r1, r2, r3 = pair_values[0]
        ra, rb, rc, rd = pair_values[1]
        # Get slope, intercept from 0-2 edge line and 1-3 edge line

        def get_line_equation(rA, rB):
            xA, yA, zA = rA
            xB, yB, zB = rB
            m = (yB - yA) / (xB - xA)
            c = yA - m * xA
            return (m, c)

        def get_line_intersection(m1, c1, m2, c2):
            x = (c2 - c1) / (m1 - m2)
            y = m1 * x + c1
            return (x, y)

        m03, c03 = get_line_equation(r0, r3)
        m12, c12 = get_line_equation(r1, r2)

        # Get slope, intercept for a-b and c-d
        mab, cab = get_line_equation(ra, rb)
        mcd, ccd = get_line_equation(rc, rd)
        # Get intersection between 0-2 with a-b
        r03_ab = get_line_intersection(m03, c03, mab, cab)
        # Get line intersection between 1-3 and a-b
        r12_ab = get_line_intersection(m12, c12, mab, cab)
        # Get line intersection between 0-2 and c-d
        r03_cd = get_line_intersection(m03, c03, mcd, ccd)
        # Get line intersection between 1-3 and c-d
        r12_cd = get_line_intersection(m12, c12, mcd, ccd)
        # returns pairs of x,y
        return [
            (r03_ab[0], r03_ab[1], r12_ab[0], r12_ab[1]),
            (r03_cd[0], r03_cd[1], r12_cd[0], r12_cd[1]),
        ]


class image_kapton_correction(object):
    def __init__(
        self,
        panel_size_px=None,  #
        pixel_size_mm=None,  #
        detector_dist_mm=None,  #
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
            ), "Kapton param sigmas must be non-negative"
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
            # instantiate Kapton absorption class here
            absorption = KaptonTape_2019(
                params_version[0],
                params_version[1],
                params_version[2],
                params_version[3],
                self.wavelength_ang,
            )
            # *map(float, self.panel_size_px)) #
            detector = self.expt.detector

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
                    f_absolute, s_absolute, z_absolute = (
                        shoebox.coords().select(foreground).parts()
                    )
                    panel_number = self.reflections_sele[iref]["panel"]
                    lab_coords = detector[panel_number].get_lab_coord(
                        detector[panel_number].pixel_to_millimeter(
                            flex.vec2_double(f_absolute, s_absolute)
                        )
                    )
                    s1 = lab_coords.each_normalize()
                    # Real step right here
                    kapton_correction_vector.extend(absorption.abs_correction_flex(s1))
                    #
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
                plt.hist(data, 20)
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
            panels = [p for p in detector]
            panel_size_px = [p.get_image_size() for p in panels]
            pixel_size_mm = [p.get_pixel_size()[0] for p in panels]
            detector_dist_mm = [p.get_distance() for p in panels]
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
