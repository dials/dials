from __future__ import annotations

import copy

import libtbx.phil
from scitbx import matrix

from dials.array_family import flex

phil_scope = libtbx.phil.parse(
    """
reverse_phi = False
  .type = bool
  .optional = True
crystal_frame = False
  .type = bool
  .optional = True
beam_centre_panel = None
  .type = int
  .help = "Panel number for the beam centre"
beam_centre = None
  .type = floats(size=2)
  .help = "Fast, slow beam centre coordinates (mm)."
d_min = None
  .type = float(value_min=0.0)
z_min = None
  .type = float
z_max = None
  .type = float
n_min = 0
  .type = int
  .help = "Minimum size of spot in pixels"
n_max = 1000000
  .type = int
  .help = "Maximum size of spot in pixels"
partiality_min = None
  .type = float
partiality_max = None
  .type = float
display = *all unindexed indexed integrated
  .type = choice
outlier_display = outliers inliers
  .type = choice
experiment_ids = None
  .type = ints(value_min=-1)
black_background = True
  .type = bool
  .help = "Switch between black or white background"
"""
)


class Render3d:
    def __init__(self, settings=None):
        self.reflections = None
        if settings is None:
            self.settings = phil_scope.fetch().extract()
        else:
            self.settings = settings
        self.goniometer_orig = None
        self.viewer = None

    def load_models(self, experiments, reflections):
        self.experiments = experiments
        self.reflections_input = reflections
        if self.experiments[0].goniometer is not None:
            self.viewer.set_rotation_axis(
                self.experiments[0].goniometer.get_rotation_axis()
            )
        self.viewer.set_beam_vector(self.experiments[0].beam.get_unit_s0())
        if self.settings.beam_centre is None:
            try:
                (
                    self.settings.beam_centre_panel,
                    self.settings.beam_centre,
                ) = self.experiments[0].detector.get_ray_intersection(
                    self.experiments[0].beam.get_s0()
                )
            except RuntimeError:
                pass
        else:
            self.set_beam_centre(
                self.settings.beam_centre_panel, self.settings.beam_centre
            )
        crystals = [
            expt.crystal for expt in self.experiments if expt.crystal is not None
        ]
        if crystals:
            # the points are scaled by 100 so must do that here too
            vecs = [
                [
                    matrix.col(l)
                    for l in (100 * matrix.sqr(c.get_A()))
                    .transpose()
                    .as_list_of_lists()
                ]
                for c in crystals
            ]
            self.viewer.set_reciprocal_lattice_vectors(vecs)
            vecs = [
                [
                    matrix.col(l)
                    for l in (100 * matrix.sqr(c.get_B()))
                    .transpose()
                    .as_list_of_lists()
                ]
                for c in crystals
            ]
            self.viewer.set_reciprocal_crystal_vectors(vecs)
        self.map_points_to_reciprocal_space()
        self.set_points()

    def map_points_to_reciprocal_space(self):
        # 155 handle data from predictions *only* if that is what we have
        calculated = "xyzobs.px.value" not in self.reflections_input
        self.reflections = copy.deepcopy(self.reflections_input)
        experiments = copy.deepcopy(self.experiments)

        if not calculated:
            self.reflections.centroid_px_to_mm(experiments)

        if self.settings.reverse_phi:
            for expt in experiments:
                expt.goniometer.set_rotation_axis(
                    [-i for i in expt.goniometer.get_rotation_axis()]
                )

        self.reflections.map_centroids_to_reciprocal_space(
            experiments,
            calculated=calculated,
            crystal_frame=self.settings.crystal_frame,
        )

    def set_points(self):
        reflections = self.reflections

        if "miller_index" in reflections:
            if "flags" not in reflections:
                reflections.set_flags(
                    reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
                )
                reflections["id"].set_selected(
                    reflections["miller_index"] == (0, 0, 0), -1
                )
                reflections.set_flags(
                    flex.bool(len(reflections), True), reflections.flags.strong
                )
                reflections.set_flags(
                    flex.bool(len(reflections), False), reflections.flags.integrated
                )
                reflections.set_flags(
                    flex.bool(len(reflections), False),
                    reflections.flags.centroid_outlier,
                )

            outlier_sel = reflections.get_flags(reflections.flags.centroid_outlier)

            if self.settings.outlier_display == "outliers":
                reflections = reflections.select(outlier_sel)
            if self.settings.outlier_display == "inliers":
                reflections = reflections.select(~outlier_sel)

            indexed_sel = reflections.get_flags(reflections.flags.indexed)
            strong_sel = reflections.get_flags(reflections.flags.strong)
            integrated_sel = reflections.get_flags(reflections.flags.integrated)

            if self.settings.display == "indexed":
                reflections = reflections.select(indexed_sel)
            elif self.settings.display == "unindexed":
                reflections = reflections.select(strong_sel & ~indexed_sel)
            elif self.settings.display == "integrated":
                reflections = reflections.select(integrated_sel)

        if self.settings.experiment_ids:
            sel = flex.bool(len(reflections), False)
            for i in self.settings.experiment_ids:
                sel.set_selected(reflections["id"] == i, True)
            reflections = reflections.select(sel)

        d_spacings = 1 / reflections["rlp"].norms()

        # 155 handle data from predictions *only* if that is what we have
        if "xyzobs.px.value" in self.reflections_input:
            use_column = "xyzobs.px.value"
        else:
            use_column = "xyzcal.px"

        if self.settings.d_min is not None:
            reflections = reflections.select(d_spacings >= self.settings.d_min)
        else:
            self.settings.d_min = flex.min(d_spacings)
        if self.settings.z_min is not None:
            z = reflections[use_column].parts()[2]
            reflections = reflections.select(z >= self.settings.z_min)
        else:
            z = reflections[use_column].parts()[2]
            self.settings.z_min = flex.min(z)
        if self.settings.z_max is not None:
            z = reflections[use_column].parts()[2]
            reflections = reflections.select(z <= self.settings.z_max)
        else:
            z = reflections[use_column].parts()[2]
            self.settings.z_max = flex.max(z)

        if "n_signal" in reflections:
            if self.settings.n_min is not None:
                _ns = reflections["n_signal"]
                reflections = reflections.select(_ns >= self.settings.n_min)
            else:
                _ns = reflections["n_signal"]
                self.settings.n_min = int(flex.min(_ns))

            if self.settings.n_max is not None:
                _ns = reflections["n_signal"]
                reflections = reflections.select(_ns <= self.settings.n_max)
            else:
                _ns = reflections["n_signal"]
                self.settings.n_max = int(flex.max(_ns))

        if "partiality" in reflections:
            p = reflections["partiality"]
            if self.settings.partiality_min is not None:
                sel = p >= self.settings.partiality_min
                reflections = reflections.select(sel)
                p = p.select(sel)
            else:
                self.settings.partiality_min = flex.min(p)
            if self.settings.partiality_max is not None:
                reflections = reflections.select(p <= self.settings.partiality_max)
            else:
                self.settings.partiality_max = flex.max(p)
        points = reflections["rlp"] * 100
        self.viewer.set_points(points)
        self.viewer.set_points_data(reflections)
        colors = flex.vec3_double(len(points), (1, 1, 1))

        if len(points):
            # suggested colorblind color palette
            # sorry if you have > 8 lattices!
            palette = flex.vec3_double(
                (
                    (255, 255, 255),
                    (230, 159, 0),
                    (86, 180, 233),
                    (0, 158, 115),
                    (240, 228, 66),
                    (0, 114, 178),
                    (213, 94, 0),
                    (204, 121, 167),
                )
            )
            palette *= 1 / 255
            if not self.settings.black_background:
                bkg = flex.vec3_double()
                for j in range(len(palette)):
                    bkg.append((1, 1, 1))
                palette = bkg - palette
            self.viewer.set_palette(palette)
            n = palette.size() - 1
            if (
                self.reflections.get_flags(self.reflections.flags.indexed).count(True)
                == 0
            ):
                if "imageset_id" in reflections:
                    imageset_id = reflections["imageset_id"]
                else:
                    imageset_id = reflections["id"]
                for i in range(0, flex.max(imageset_id) + 1):
                    colors.set_selected(imageset_id == i, palette[(i % n) + 1])
            else:
                colors.set_selected(reflections["id"] == -1, palette[0])
                for i in range(0, flex.max(reflections["id"]) + 1):
                    colors.set_selected(reflections["id"] == i, palette[(i % n) + 1])
        self.viewer.set_colors(colors)

    def set_beam_centre(self, panel, beam_centre):
        detector = self.experiments[0].detector
        beam = self.experiments[0].beam
        panel = int(panel)

        try:
            p = detector[panel]
        except RuntimeError:
            raise ValueError(f"Detector does not have panel index {panel}")

        # Check if the new beam centre can be set
        beam_f, beam_s = beam_centre
        trial_s0_direction = p.get_lab_coord((beam_f, beam_s))
        try:
            panel_id, beam_centre = detector.get_ray_intersection(trial_s0_direction)
        except RuntimeError:
            # beam centre calculation fails if the beam falls between panels
            raise ValueError("No detector intersection")

        if panel_id != panel:
            raise ValueError(f"Beam centre cannot be set with panel {panel}")

        # If we get this far, then the beam centre is valid for this panel
        beam.set_unit_s0(trial_s0_direction)
