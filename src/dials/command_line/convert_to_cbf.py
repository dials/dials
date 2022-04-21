from __future__ import annotations

import sys

import iotbx.phil

import dials.util

help_message = """

Convert data which can be read by DIALS, given a experiment list, to CBF format -
with ersatz miniCBF header. Can be used e.g. with HDF5 format data.

Examples::

  dials.convert_to_cbf models.expt prefix=data_as_cbf
"""

phil_scope = iotbx.phil.parse(
    """\
output {
  template = as_cbf_%04d.cbf
    .type = path
}
""",
    process_includes=True,
)


def convert_to_cbf(imageset, template):
    from dxtbx.format.FormatCBFMini import FormatCBFMini

    for i in range(len(imageset)):
        print(template % (i + 1))
        FormatCBFMini.as_file(
            imageset.get_detector(),
            imageset.get_beam(),
            imageset.get_goniometer(),
            imageset.get_scan()[i],
            imageset.get_raw_data(i)[0],
            template % (i + 1),
        )


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util.options import ArgumentParser, flatten_experiments

    usage = "dials.convert_to_cbf [options] models.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args, show_diff_phil=True, return_unhandled=True
    )

    template = params.output.template
    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        return

    if len(experiments) > 1:
        sys.exit("Only one experiment can be processed at a time")
    else:
        imagesets = experiments.imagesets()
        assert len(imagesets) == 1, len(imagesets)
        imageset = imagesets[0]

    convert_to_cbf(imageset, template)


if __name__ == "__main__":
    run()
