import iotbx.phil

import dials.util
from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

help_message = """
Small modifications to an indexed spot list to allow it to be treated as if
it were integrated, for rapid feedback analysis.
"""

phil_scope = iotbx.phil.parse(
    """
output {
  reflections = indexed_plus.refl
    .type = path
}
"""
)


def indexed_as_integrated(reflections, params, experiments):
    """Small modifications to an indexed spot list to allow it to be
    treated as if it were integrated, for rapid feedback analysis."""

    # filter only indexed reflections & assign the summation integrated flag

    sel = reflections.get_flags(reflections.flags.indexed)
    reflections = reflections.select(sel)
    all = flex.bool(reflections.size(), True)
    reflections.set_flags(all, reflections.flags.integrated_sum)

    # add resolution to reflections

    if "imageset_id" not in reflections:
        reflections["imageset_id"] = reflections["id"]

    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    reflections["d"] = 1 / reflections["rlp"].norms()

    return reflections


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.indexed_as_integrated [options] indexed.refl indexed.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(reflections) != 1 or len(experiments) != 1:
        parser.print_help()
        return

    pseudo_integrated = indexed_as_integrated(reflections[0], params, experiments)
    pseudo_integrated.as_file(params.output.reflections)


if __name__ == "__main__":
    run()
