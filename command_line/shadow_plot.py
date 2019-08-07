# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import absolute_import, division, print_function

import json
import sys

import libtbx
import libtbx.phil
from dials.util import Sorry
from scitbx.array_family import flex

help_message = """
Generate a 1d or 2d goniometer detector shadow plot for a given experiment list.

Examples::

  dials.shadow_plot models.expt

  dials.shadow_plot models.expt mode=2d

"""

phil_scope = libtbx.phil.parse(
    """
oscillation_range = None
  .type = floats(size=2)
step_size = auto
  .type = float(value_min=0)
y_max = None
  .type = float(value_min=0, value_max=100)
mode = *1d 2d
  .type = choice
output {
  plot = scan_shadow_plot.png
    .type = path
  json = None
    .type = path
  size_inches = None
    .type = floats(value_min=0, size=2)
}
"""
)


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments

    usage = "dials.shadow_plot [options] models.expt"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        sys.exit(0)

    assert len(experiments) == 1
    imagesets = experiments.imagesets()

    imageset = imagesets[0]
    goniometer = imageset.get_goniometer()
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    masker = imageset.masker()
    if masker is None:
        raise Sorry("Goniometer model does not support shadowing.")
    angles = goniometer.get_angles()
    names = goniometer.get_names()
    scan_axis = goniometer.get_scan_axis()
    phi = angles[0]

    if params.step_size is libtbx.Auto:
        if params.mode == "1d":
            step = scan.get_oscillation()[1]
        else:
            step = 10
    else:
        step = params.step_size

    if params.mode == "1d":
        if params.oscillation_range is not None:
            start, end = params.oscillation_range
        else:
            start, end = scan.get_oscillation_range()

        scan_points = flex.double(libtbx.utils.frange(start, end, step=step))
        n_px_shadowed = flex.double(scan_points.size(), 0)
        n_px_tot = flex.double(scan_points.size(), 0)

        assert len(angles) == 3
        for i, scan_angle in enumerate(scan_points):
            shadow = masker.project_extrema(detector, scan_angle)
            for p_id in range(len(detector)):
                px_x, px_y = detector[p_id].get_image_size()
                n_px_tot[i] += px_x * px_y
                if shadow[p_id].size() < 4:
                    continue
                n_px_shadowed[i] += polygon_area(shadow[p_id])

    else:
        kappa_values = flex.double(libtbx.utils.frange(0, 360, step=step))
        omega_values = flex.double(libtbx.utils.frange(0, 360, step=step))
        grid = flex.grid(kappa_values.size(), omega_values.size())
        n_px_shadowed = flex.double(grid, 0)
        n_px_tot = flex.double(grid, 0)

        assert len(angles) == 3
        for i, kappa in enumerate(kappa_values):
            for j, omega in enumerate(omega_values):
                masker.set_goniometer_angles((phi, kappa, omega))
                masker.extrema_at_scan_angle(omega)
                shadow = masker.project_extrema(detector, omega)
                for p_id in range(len(detector)):
                    px_x, px_y = detector[p_id].get_image_size()
                    n_px_tot[i, j] += px_x * px_y
                    if shadow[p_id].size() < 4:
                        continue
                    n_px_shadowed[i, j] += polygon_area(shadow[p_id])

    fraction_shadowed = n_px_shadowed / n_px_tot

    if params.output.json is not None:
        if params.mode == "2d":
            raise Sorry("json output not supported for mode=2d")

        print("Writing json output to %s" % params.output.json)
        d = {
            "scan_points": list(scan_points),
            "fraction_shadowed": list(fraction_shadowed),
        }
        with open(params.output.json, "w") as f:
            json.dump(d, f)

    if params.output.plot is not None:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt

        plt.style.use("ggplot")

        if params.mode == "1d":
            plt.plot(
                scan_points.as_numpy_array(), fraction_shadowed.as_numpy_array() * 100
            )
            plt.xlabel("%s angle (degrees)" % names[scan_axis])
            plt.ylabel("Shadowed area (%)")
            if params.y_max is not None:
                plt.ylim(0, params.y_max)
            else:
                plt.ylim(0, plt.ylim()[1])
        else:
            fig = plt.imshow(
                fraction_shadowed.as_numpy_array() * 100, interpolation="bicubic"
            )
            plt.xlabel("%s angle (degrees)" % names[2])
            plt.ylabel("%s angle (degrees)" % names[1])
            plt.xlim(0, 360 / step)
            plt.ylim(0, 360 / step)
            fig.axes.set_xticklabels(["%.0f" % (step * t) for t in plt.xticks()[0]])
            fig.axes.set_yticklabels(["%.0f" % (step * t) for t in plt.yticks()[0]])
            cbar = plt.colorbar()
            cbar.set_label("Shadowed area (%)")

        if params.output.size_inches is not None:
            fig = plt.gcf()
            fig.set_size_inches(params.output.size_inches)
        plt.tight_layout()
        print("Saving plot to %s" % params.output.plot)
        plt.savefig(params.output.plot)


def polygon_area(points):
    # http://mathworld.wolfram.com/PolygonArea.html
    x0, y0 = points.parts()
    x1 = x0[1:]
    x1.append(x0[0])
    y1 = y0[1:]
    y1.append(y0[0])

    return 0.5 * abs(flex.sum(x0 * y1 - x1 * y0))


if __name__ == "__main__":
    run(sys.argv[1:])
