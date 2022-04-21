# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from __future__ import annotations

import logging
import math

import libtbx.phil
from scitbx import matrix
from scitbx.array_family import flex

import dials.util
from dials.algorithms.indexing.indexer import find_max_cell
from dials.command_line.search_beam_position import run_dps
from dials.util.reciprocal_lattice import Render3d

try:
    import matplotlib
except ImportError:
    exit()  # To pass through the "make" step, for graphics-free HPC build

# Offline backend
matplotlib.use("Agg")

logger = logging.getLogger("dials.command_line.rl_png")

help_message = """
Generate a png of the strong spots from spotfinding in reciprocal space.

Examples::

  dials.rl_png imported.expt strong.refl

  dials.rl_png indexed.expt indexed.refl
"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.reciprocal_lattice.phil_scope

marker_size = 5
  .type = int(value_min=1)

basis_vector_search {
  n_solutions = 3
    .type = int
}
plot {
  size_inches = 10,10
    .type = floats(size=2, value_min=0)
}
""",
    process_includes=True,
)


class ReciprocalLatticePng(Render3d):
    def __init__(self, settings=None):
        Render3d.__init__(self, settings=settings)
        self.viewer = PngScene(settings=self.settings)


class PngScene:
    def __init__(self, settings):
        self.settings = settings
        self.rotation_axis = None
        self.beam_vector = None
        self.points = None
        self.colors = None
        self.palette = None

    def set_rotation_axis(self, axis):
        self.rotation_axis = axis

    def set_beam_vector(self, beam):
        self.beam_vector = beam

    def set_points(self, points):
        self.points = points

    def set_points_data(self, reflections):
        # we do not label reciprocal lattice points here
        pass

    def set_colors(self, colors):
        # convert whites to black (background is white)
        colors.set_selected((colors.norms() == math.sqrt(3)), (0, 0, 0))
        self.colors = colors

    def set_palette(self, palette):
        self.palette = palette

    def set_reciprocal_lattice_vectors(self, *args, **kwargs):
        # we do not draw reciprocal lattice vectors at this time
        pass

    def set_reciprocal_crystal_vectors(self, *args, **kwargs):
        # we do not draw reciprocal crystal vectors at this time either
        pass

    def project_2d(self, n):
        d = self.points.dot(n.elems)
        p = d * flex.vec3_double(len(d), n.elems)

        points2d = self.points - p

        x = matrix.col((1, 0, 0))
        if x.angle(n) == 0 or x.angle(-n) == 0:
            x = matrix.col((0, 1, 0))

        x = (x - x.dot(n) * n).normalize()

        y = x.cross(n)
        # assert y.angle(x, deg=True) == 90
        # assert y.angle(matrix.col(n), deg=True) == 90

        px2d = points2d.dot(x)
        py2d = points2d.dot(y)

        return px2d, py2d

    def plot(self, filename, n=(1, 0, 0)):
        from matplotlib import pyplot

        n = matrix.col(n).normalize()
        x, y = self.project_2d(n)
        fig = pyplot.figure(figsize=self.settings.plot.size_inches)
        pyplot.scatter(
            x.as_numpy_array(),
            y.as_numpy_array(),
            marker="+",
            s=self.settings.marker_size,
            c=list(self.colors),
        )
        pyplot.title("Plane normal: (%.2g, %.2g, %.2g)" % (n.elems))
        fig.savefig(filename)
        pyplot.close()


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util import log
    from dials.util.options import (
        ArgumentParser,
        reflections_and_experiments_from_files,
    )

    usage = "dials.rl_png [options] experiments.json observations.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    # Configure the logging
    log.config(logfile="dials.rl_png.log")

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections = reflections[0]

    f = ReciprocalLatticePng(settings=params)
    f.load_models(experiments, reflections)

    rotation_axis = matrix.col(experiments[0].goniometer.get_rotation_axis())
    s0 = matrix.col(experiments[0].beam.get_s0())

    e1 = rotation_axis.normalize()
    e2 = s0.normalize()
    e3 = e1.cross(e2).normalize()

    f.viewer.plot("rl_rotation_axis.png", n=e1.elems)
    f.viewer.plot("rl_beam_vector", n=e2.elems)
    f.viewer.plot("rl_e3.png", n=e3.elems)

    n_solutions = params.basis_vector_search.n_solutions

    if experiments.crystals().count(None) < len(experiments):
        for i, c in enumerate(experiments.crystals()):
            A = matrix.sqr(c.get_A())

            direct_matrix = A.inverse()
            a = direct_matrix[:3]
            b = direct_matrix[3:6]
            c = direct_matrix[6:9]

            prefix = ""
            if len(experiments.crystals()) > 1:
                prefix = "%i_" % (i + 1)

            f.viewer.plot(f"rl_{prefix}a.png", n=a)
            f.viewer.plot(f"rl_{prefix}b.png", n=b)
            f.viewer.plot(f"rl_{prefix}c.png", n=c)

    elif n_solutions:
        if "imageset_id" not in reflections:
            reflections["imageset_id"] = reflections["id"]

        reflections.centroid_px_to_mm(experiments)

        reflections.map_centroids_to_reciprocal_space(experiments)

        if params.d_min is not None:
            d_spacings = 1 / reflections["rlp"].norms()
            sel = d_spacings > params.d_min
            reflections = reflections.select(sel)

        # derive a max_cell from mm spots
        max_cell = find_max_cell(
            reflections, max_cell_multiplier=1.3, step_size=45
        ).max_cell

        result = run_dps(experiments[0], reflections, max_cell)
        if result:
            solutions = [matrix.col(v) for v in result["solutions"]]
            for i in range(min(n_solutions, len(solutions))):
                v = solutions[i]
                f.viewer.plot(f"rl_solution_{i + 1}.png", n=v.elems)


if __name__ == "__main__":
    run()
