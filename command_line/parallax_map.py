from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME dev.dials.parallax_map

import libtbx.phil
from scitbx.array_family import flex
import math

scope = libtbx.phil.parse(
    """
  output = map.png
    .type = path
"""
)


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    import libtbx.load_env

    usage = "%s [options] experiments.json" % (libtbx.env.dispatcher_name)

    parser = OptionParser(usage=usage, phil=scope, read_experiments=True)

    params, options = parser.parse_args(show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) != 1:
        raise RuntimeError("please pass exactly one experiment")

    experiment = experiments[0]

    panel = experiment.detector[0]

    # we know the nominal pixel size, so compute a map of how far the shift
    # between the apparent and effective pixel position is

    ny, nx = panel.get_image_size()
    dy, dx = panel.get_pixel_size()

    d_map = flex.double(flex.grid((ny, nx)))

    for j in range(ny):
        for i in range(nx):
            _x, _y = panel.pixel_to_millimeter((i, j))
            _i = _x / dx
            _j = _y / dy
            d_map[j, i] = math.sqrt((i - _i) ** 2 + (j - _j) ** 2)

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot

    pyplot.imshow(d_map.as_numpy_array())
    pyplot.colorbar()
    pyplot.savefig(params.output)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
