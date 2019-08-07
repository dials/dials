from __future__ import absolute_import, division, print_function

import libtbx.phil

help_message = """

"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.options.geometry_phil_scope
output {
  experiments = modified.expt
    .type = path
}
""",
    process_includes=True,
)


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments

    usage = "dials.modify_geometry [options] models.expt"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    from dials.command_line.dials_import import ManualGeometryUpdater

    update_geometry = ManualGeometryUpdater(params)

    if len(experiments):
        imagesets = experiments.imagesets()

    for imageset in imagesets:
        imageset_new = update_geometry(imageset)
        imageset.set_detector(imageset_new.get_detector())
        imageset.set_beam(imageset_new.get_beam())
        imageset.set_goniometer(imageset_new.get_goniometer())
        imageset.set_scan(imageset_new.get_scan())

    from dxtbx.serialize import dump

    if len(experiments):
        print("Saving modified experiments to %s" % params.output.experiments)
        dump.experiment_list(experiments, params.output.experiments)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
