# LIBTBX_SET_DISPATCHER_NAME dials.modify_geometry
# LIBTBX_SET_DISPATCHER_NAME dials.modify_experiment

from __future__ import absolute_import, division, print_function

import warnings

import libtbx.phil

import dials.util

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


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util.options import OptionParser, flatten_experiments

    usage = "dials.modify_geometry [options] models.expt"

    import libtbx.load_env

    if libtbx.env.dispatcher_name == "dials.modify_geometry":
        warnings.warn(
            "dials.modify_geometry is now deprecated, please use dials.modify_experiment instead"
        )

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)
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

    if len(experiments):
        print("Saving modified experiments to %s" % params.output.experiments)
        experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()
