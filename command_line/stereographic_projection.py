# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
from __future__ import absolute_import, division, print_function

import json
import math
import os
import sys

import iotbx.phil
from cctbx import crystal, miller
from cctbx.array_family import flex
from scitbx import matrix

help_message = """

Calculates a stereographic projection image for the given crystal models and
the given miller indices (either specified invidually, or for all miller indices
up to a given hkl_limit). By default the projection is in the plane
perpendicular to 0,0,1 reflection for the first crystal, however the projection
can optionally be performed in the laboratory frame (frame=laboratory) in the
plane perpendicular to the beam. Setting the parameter expand_to_p1=True will
also plot all symmetry equivalents of the given miller indices, and
eliminate_sys_absent=False will eliminate systematically absent reflections
before generating the projection.

Examples::

  dials.stereographic_projection experiments.json hkl=1,0,0 hkl=0,1,0

  dials.stereographic_projection experiments.json hkl_limit=2

  dials.stereographic_projection experiments_1.json experiments_2.json hkl=1,0,0 expand_to_p1=True

"""

phil_scope = iotbx.phil.parse(
    """
hkl = None
  .type = ints(size=3)
  .multiple=True
hkl_limit = None
  .type = int(value_min=1)
expand_to_p1 = True
  .type = bool
  .help = "Expand the given miller indices to symmetry equivalent reflections"
eliminate_sys_absent = False
  .type = bool
  .help = "Eliminate systematically absent reflections"
frame = *laboratory crystal
  .type = choice
phi_angle = 0
  .type = float
  .help = "Phi rotation angle (degrees)"
use_starting_angle = False
  .type = bool
  .help = "If True, then the projection will be done for each crystal at the "
          "starting phi angle for the scan associated with the crystal."
plane_normal = None
  .type = ints(size=3)
save_coordinates = True
  .type = bool
plot {
  show = False
    .type = bool
  filename = stereographic_projection.png
    .type = path
  label_indices = False
    .type = bool
  colours = None
    .type = strings
  marker_size = 3
    .type = int(value_min=1)
  font_size = 6
    .type = float(value_min=0)
  colour_map = None
    .type = str
  gridsize = None
    .type = int
}
json {
  filename = None
    .type = path
}
"""
)


def reference_poles_perpendicular_to_beam(beam, goniometer):
    # plane normal
    d0 = matrix.col(beam.get_s0()).normalize()
    if goniometer is not None:
        d1 = d0.cross(matrix.col(goniometer.get_rotation_axis())).normalize()
    else:
        d1 = d0.ortho()
    d2 = d1.cross(d0).normalize()
    return (d0, d1, d2)


def reference_poles_crystal(crystal_model, plane_normal=(0, 0, 1)):
    A = matrix.sqr(crystal_model.get_A())
    B = matrix.sqr(crystal_model.get_B())
    A_inv = A.inverse()
    G = A_inv * A_inv.transpose()
    G_star = A.transpose() * A
    h0 = (G * matrix.col(plane_normal)).normalize()
    h1 = matrix.col((1, 0, 0)).cross((G_star * h0).normalize())
    h2 = (G_star * h1).cross(G_star * h0).normalize()
    return tuple((B * h).normalize() for h in (h0, h1, h2))


def stereographic_projection(points, reference_poles):
    # https://doi.org/10.1107/S0021889868005029
    # J. Appl. Cryst. (1968). 1, 68-70
    # The construction of stereographic projections by computer
    # G. K. Stokes, S. R. Keown and D. J. Dyson

    assert len(reference_poles) == 3
    r_0, r_1, r_2 = reference_poles

    projections = flex.vec2_double()

    for p in points:
        r_i = matrix.col(p)
        # theta is the angle between r_i and the plane normal, r_0
        cos_theta = r_i.cos_angle(r_0)
        if cos_theta < 0:
            r_i = -r_i
            cos_theta = r_i.cos_angle(r_0)

        # alpha is the angle between r_i and r_1
        cos_alpha = r_i.cos_angle(r_1)
        theta = math.acos(cos_theta)
        cos_phi = cos_alpha / math.sin(theta)
        if abs(cos_phi) > 1:
            cos_phi = math.copysign(1, cos_phi)
        phi = math.acos(cos_phi)

        N = r_i.dot(r_2)
        r = math.tan(theta / 2)
        x = r * cos_phi
        y = r * math.sin(phi)
        y = math.copysign(y, N)

        projections.append((x, y))

    return projections


def gcd_list(l):
    # greatest common divisor for a list of numbers
    from scitbx.math import gcd_int_simple as gcd

    result = l[0]
    for i in range(1, len(l)):
        result = gcd(result, l[i])
    return result


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments

    # The script usage
    usage = "dials.stereographic_projection [options] [param.phil] experiments.json"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    if not experiments:
        parser.print_help()
        return

    if not params.hkl and params.hkl_limit is None:
        sys.exit("Please provide hkl or hkl_limit parameters.")

    if params.hkl is not None and len(params.hkl):
        miller_indices = flex.miller_index(params.hkl)
    elif params.hkl_limit is not None:
        limit = params.hkl_limit
        miller_indices = flex.miller_index()
        for h in range(-limit, limit + 1):
            for k in range(-limit, limit + 1):
                for l in range(-limit, limit + 1):
                    if (h, k, l) == (0, 0, 0):
                        continue
                    miller_indices.append((h, k, l))

    crystals = experiments.crystals()

    symmetry = crystal.symmetry(
        unit_cell=crystals[0].get_unit_cell(), space_group=crystals[0].get_space_group()
    )
    miller_set = miller.set(symmetry, miller_indices)
    d_spacings = miller_set.d_spacings()
    if params.eliminate_sys_absent:
        d_spacings = d_spacings.eliminate_sys_absent()
    if params.expand_to_p1:
        d_spacings = d_spacings.as_non_anomalous_array().expand_to_p1()
        d_spacings = d_spacings.generate_bijvoet_mates()
    miller_indices = d_spacings.indices()

    # find the greatest common factor (divisor) between miller indices
    miller_indices_unique = flex.miller_index()
    for hkl in miller_indices:
        gcd = gcd_list(hkl)
        if gcd > 1:
            miller_indices_unique.append(tuple(int(h / gcd) for h in hkl))
        elif gcd < 1:
            pass
        else:
            miller_indices_unique.append(hkl)
    miller_indices = miller_indices_unique
    miller_indices = flex.miller_index(list(set(miller_indices)))

    ref_crystal = crystals[0]
    U = matrix.sqr(ref_crystal.get_U())
    B = matrix.sqr(ref_crystal.get_B())
    R = matrix.identity(3)

    if params.frame == "laboratory":
        reference_poles = reference_poles_perpendicular_to_beam(
            experiments[0].beam, experiments[0].goniometer
        )
        if params.use_starting_angle:
            rotation_axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
            R = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
                experiments[0].scan.get_oscillation()[0], deg=True
            )
        elif params.phi_angle != 0:
            rotation_axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
            R = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
                params.phi_angle, deg=True
            )
    else:
        if params.plane_normal is not None:
            plane_normal = params.plane_normal
        else:
            plane_normal = (0, 0, 1)
        reference_poles = reference_poles_crystal(
            ref_crystal, plane_normal=plane_normal
        )

    if params.frame == "crystal":
        U = matrix.identity(3)

    reciprocal_space_points = list(R * U * B) * miller_indices.as_vec3_double()
    projections_ref = stereographic_projection(reciprocal_space_points, reference_poles)

    projections_all = [projections_ref]

    if experiments:
        from dials.algorithms.indexing.compare_orientation_matrices import (
            difference_rotation_matrix_axis_angle,
        )

        for expt in experiments[1:]:
            cryst = expt.crystal
            if params.frame == "crystal":
                R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
                    ref_crystal, cryst
                )
                U = R_ij
            elif params.use_starting_angle:
                if params.use_starting_angle:
                    rotation_axis = matrix.col(expt.goniometer.get_rotation_axis())
                    R = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
                        expt.scan.get_oscillation()[0], deg=True
                    )
            else:
                U = matrix.sqr(cryst.get_U())
            reciprocal_space_points = (
                list(R * U * matrix.sqr(cryst.get_B()))
                * miller_indices.as_vec3_double()
            )
            projections = stereographic_projection(
                reciprocal_space_points, reference_poles
            )
            projections_all.append(projections)

    if params.save_coordinates:
        with open("projections.txt", "w") as f:
            f.write("crystal h k l x y" + os.linesep)
            for i_cryst, projections in enumerate(projections_all):
                for hkl, proj in zip(miller_indices, projections):
                    f.write("%i " % (i_cryst + 1))
                    f.write("%i %i %i " % hkl)
                    f.write(("%f %f" + os.linesep) % proj)

    if params.plot.show or params.plot.filename:
        epochs = None
        if params.plot.colour_map is not None:
            if experiments[0].scan is not None:
                epochs = [expt.scan.get_epochs()[0] for expt in experiments]
            else:
                epochs = [i for i, expt in enumerate(experiments)]
        plot_projections(
            projections_all,
            filename=params.plot.filename,
            show=params.plot.show,
            colours=params.plot.colours,
            marker_size=params.plot.marker_size,
            font_size=params.plot.font_size,
            gridsize=params.plot.gridsize,
            label_indices=params.plot.label_indices,
            epochs=epochs,
            colour_map=params.plot.colour_map,
        )

    if params.json.filename:
        projections_as_json(projections_all, filename=params.json.filename)


def plot_projections(
    projections,
    filename=None,
    show=None,
    colours=None,
    marker_size=3,
    font_size=6,
    gridsize=None,
    label_indices=False,
    epochs=None,
    colour_map=None,
):
    assert [filename, show].count(None) < 2
    projections_all = projections

    try:
        import matplotlib

        if not show:
            # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
            matplotlib.use("Agg")  # use a non-interactive backend
        from matplotlib import pyplot
        from matplotlib import pylab
    except ImportError:
        raise Sorry("matplotlib must be installed to generate a plot.")

    if epochs is not None and colour_map is not None:
        epochs = flex.double(epochs)
        epochs -= flex.min(epochs)
        epochs /= flex.max(epochs)
        cmap = matplotlib.cm.get_cmap(colour_map)
        colours = [cmap(e) for e in epochs]
    elif colours is None or len(colours) == 0:
        colours = ["b"] * len(projections_all)
    elif len(colours) < len(projections_all):
        colours = colours * len(projections_all)

    fig = pyplot.figure()

    pyplot.scatter([0], [0], marker="+", c="0.75", s=100)
    cir = pylab.Circle((0, 0), radius=1.0, fill=False, color="0.75")
    pylab.gca().add_patch(cir)

    if gridsize is not None:
        x = flex.double()
        y = flex.double()
        for i, projections in enumerate(projections_all):
            x_, y_ = projections.parts()
            x.extend(x_)
            y.extend(y_)
        hb = pyplot.hexbin(x, y, gridsize=gridsize, linewidths=0.2)
        pyplot.colorbar(hb)
    else:
        for i, projections in enumerate(projections_all):
            x, y = projections.parts()
            pyplot.scatter(
                x.as_numpy_array(),
                y.as_numpy_array(),
                c=colours[i],
                s=marker_size,
                edgecolors="none",
            )
            if label_indices:
                for j, (hkl, proj) in enumerate(zip(miller_indices, projections)):
                    # hack to not write two labels on top of each other
                    p1, p2 = (projections - proj).parts()
                    if (flex.sqrt(flex.pow2(p1) + flex.pow2(p2)) < 1e-3).iselection()[
                        0
                    ] != j:
                        continue
                    pyplot.text(proj[0], proj[1], str(hkl), fontsize=font_size)
    fig.axes[0].set_aspect("equal")
    pyplot.xlim(-1.1, 1.1)
    pyplot.ylim(-1.1, 1.1)
    if filename is not None:
        pyplot.savefig(filename, size_inches=(24, 18), dpi=300)
    if show:
        pyplot.show()


def projections_as_dict(projections):
    projections_all = flex.vec2_double()
    for proj in projections:
        projections_all.extend(proj)

    data = []
    x, y = projections_all.parts()
    data.append(
        {
            "x": list(x),
            "y": list(y),
            "mode": "markers",
            "type": "scatter",
            "name": "stereographic_projections",
            "showlegend": False,
        }
    )
    data.append(
        {
            "x": [0],
            "y": [0],
            "mode": "markers",
            "marker": {
                "color": "black",
                "size": 25,
                "symbol": "cross-thin",
                "line": {"width": 1},
            },
            "showlegend": False,
        }
    )

    d = {
        "data": data,
        "layout": {
            "title": "Stereographic projections",
            "hovermode": False,
            "xaxis": {
                "range": [-1.0, 1.0],
                "showgrid": False,
                "zeroline": False,
                "showline": False,
                "ticks": "",
                "showticklabels": False,
            },
            "yaxis": {
                "range": [-1.0, 1.0],
                "showgrid": False,
                "zeroline": False,
                "showline": False,
                "ticks": "",
                "showticklabels": False,
            },
            "shapes": [
                {
                    "type": "circle",
                    "xref": "x",
                    "yref": "y",
                    "x0": -1,
                    "y0": -1,
                    "x1": 1,
                    "y1": 1,
                    "line": {"color": "black"},
                }
            ],
        },
    }
    return d


def projections_as_json(projections, filename=None):
    d = projections_as_dict(projections)

    json_str = json.dumps(d)
    if filename is not None:
        with open(filename, "w") as f:
            f.write(json_str)
    return json_str


if __name__ == "__main__":
    run(sys.argv[1:])
