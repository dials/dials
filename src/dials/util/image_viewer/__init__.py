from __future__ import annotations

import math

from cctbx import uctbx
from dxtbx.model import Beam, Detector
from iotbx.detectors import FlexImage, FlexImage_d
from scitbx import matrix

from dials.array_family import flex

__all__ = ["calculate_isoresolution_lines"]


def calculate_isoresolution_lines(
    spacings: flex.double,
    beam: Beam,
    detector: Detector,
    flex_image: FlexImage | FlexImage_d,
    add_text: bool = True,
    n_rays: int = 720,
    binning: int = 1,
) -> tuple[
    list[tuple[tuple[float, float], tuple[float, float]]],
    list[tuple[float, float, str]],
]:
    # Calculate 2θ angles
    wavelength = beam.get_wavelength()
    twotheta = uctbx.d_star_sq_as_two_theta(uctbx.d_as_d_star_sq(spacings), wavelength)

    # Get beam vector and two orthogonal vectors
    beamvec = matrix.col(beam.get_s0())
    bor1 = beamvec.ortho()
    bor2 = beamvec.cross(bor1)

    ring_data = []
    resolution_text_data = []
    for tt, d in zip(twotheta, spacings):
        # Generate rays at 2θ
        cone_base_centre = beamvec * math.cos(tt)
        cone_base_radius = (beamvec * math.sin(tt)).length()
        rad1 = bor1.normalize() * cone_base_radius
        rad2 = bor2.normalize() * cone_base_radius
        ticks = (2 * math.pi / n_rays) * flex.double_range(n_rays)
        offset1 = flex.vec3_double(n_rays, rad1) * flex.cos(ticks)
        offset2 = flex.vec3_double(n_rays, rad2) * flex.sin(ticks)
        rays = flex.vec3_double(n_rays, cone_base_centre) + offset1 + offset2

        # Duplicate the first ray to close the loop
        rays.append(rays[0])

        # Get the ray intersections. Need to set a dummy phi value
        rt = flex.reflection_table.empty_standard(n_rays + 1)
        rt["s1"] = rays
        rt["phi"] = flex.double(n_rays + 1, 0)
        from dials.algorithms.spot_prediction import ray_intersection

        intersect = ray_intersection(detector, rt)
        rt = rt.select(intersect)
        if len(rt) == 0:
            continue

        curr_panel_id = rt[0]["panel"]
        panel = detector[curr_panel_id]

        # Split the intersections into sets of vertices in separate paths
        paths = []
        vertices: list[tuple[float, float]] = []
        for ref in rt.rows():
            if ref["panel"] != curr_panel_id:
                # close off the current path and reset the vertices
                paths.append(vertices)
                vertices = []
                curr_panel_id = ref["panel"]
                panel = detector[curr_panel_id]
            x, y = panel.millimeter_to_pixel(ref["xyzcal.mm"][0:2])
            try:
                # Multi-panel case
                y, x = flex_image.tile_readout_to_picture(
                    curr_panel_id, y - 0.5, x - 0.5
                )
            except AttributeError:
                # Single panel FlexImage
                pass
            vertices.append((x / binning, y / binning))
        paths.append(vertices)

        # For each path, convert vertices to segments and add to the ring data
        segments = []
        for vertices in paths:
            for i in range(len(vertices) - 1):
                # Avoid long segments along the image edges
                dx = abs(vertices[i + 1][0] - vertices[i][0])
                dy = abs(vertices[i + 1][1] - vertices[i][1])
                if dx > 30 or dy > 30:
                    continue
                segments.append((vertices[i], vertices[i + 1]))
        ring_data.extend(segments)

        # Add labels to the iso-resolution lines
        if add_text:
            cb1 = beamvec.rotate_around_origin(axis=bor1, angle=tt)
            for angle in (45, 135, 225, 315):
                txtvec = cb1.rotate_around_origin(
                    axis=beamvec, angle=math.radians(angle)
                )
                try:
                    panel_id, txtpos = detector.get_ray_intersection(txtvec)
                except RuntimeError:
                    continue
                txtpos = detector[panel_id].millimeter_to_pixel(txtpos)
                try:
                    # Multi-panel case
                    x, y = flex_image.tile_readout_to_picture(
                        panel_id, txtpos[1], txtpos[0]
                    )[::-1]
                except AttributeError:
                    # Single panel FlexImage
                    x, y = txtpos
                resolution_text_data.append(
                    (
                        x / binning,
                        y / binning,
                        f"{d:.2f}",
                    )
                )

    return (ring_data, resolution_text_data)
