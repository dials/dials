from __future__ import absolute_import, division, print_function

import libtbx.phil

help_message = """

Compute resolution-wise distribution of spots

Examples::

  dials.spot_resolution_shells models.expt strong.refl

"""

phil_scope = libtbx.phil.parse(
    """
  shells = 100
    .type = int(value_min=1)
"""
)


def settings():
    return phil_scope.fetch().extract()


def spot_resolution_shells(imagesets, reflections, params):
    from dials.array_family import flex

    mapped_reflections = flex.reflection_table()
    for i, imageset in enumerate(imagesets):
        if "imageset_id" in reflections:
            sel = reflections["imageset_id"] == i
        else:
            sel = reflections["id"] == i
        if isinstance(reflections["id"], flex.size_t):
            reflections["id"] = reflections["id"].as_int()
        refl = reflections.select(sel)
        refl.centroid_px_to_mm(imageset.get_detector(), imageset.get_scan())
        refl.map_centroids_to_reciprocal_space(
            imageset.get_detector(), imageset.get_beam(), imageset.get_goniometer()
        )
        mapped_reflections.extend(refl)
    reflections = mapped_reflections
    two_theta_array = reflections["rlp"].norms()
    h0 = flex.weighted_histogram(two_theta_array ** 2, n_slots=params.shells)
    n = h0.slots()
    d = 1.0 / flex.sqrt(h0.slot_centers())

    for j in range(params.shells):
        print("%d %f %d" % (j, d[j], n[j]))


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    from dials.util.options import flatten_reflections

    usage = "dials.spot_resolution_shells [options] models.expt observations.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    reflections = reflections[0]

    imagesets = experiments.imagesets()

    spot_resolution_shells(imagesets, reflections, params)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
