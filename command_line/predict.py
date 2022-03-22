from __future__ import annotations

from libtbx.phil import parse

from dials.algorithms.shadowing.filter import filter_shadowed_reflections
from dials.array_family import flex
from dials.util import show_mail_handle_errors
from dials.util.command_line import Command
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """

This program takes a set of experiments and predicts the reflections. The
reflections are then saved to file.

Examples::

  dials.predict models.expt

  dials.predict models.expt force_static=True

  dials.predict models.expt d_min=2.0
"""

phil_scope = parse(
    """
  output = predicted.refl
    .type = str
    .help = "The filename for the predicted reflections"

  force_static = False
    .type = bool
    .help = "For a scan varying model, force static prediction"

  ignore_shadows = True
    .type = bool
    .help = "Ignore dynamic shadowing"

  buffer_size = 0
    .type = int
    .help = "Calculate predictions within a buffer zone of n images either"
            " side of the scan"

  d_min = None
    .type = float
    .help = "Minimum d-spacing of predicted reflections"

    include scope dials.algorithms.profile_model.factory.phil_scope
""",
    process_includes=True,
)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = (
            "usage: dials.predict [options] [param.phil] "
            "{sequence.expt | image1.file [image2.file ...]}"
        )

        # Create the parser
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            check_format=True,
            read_experiments=True,
        )

    def run(self, args=None):
        """Execute the script."""
        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=True)

        # Check the number of experiments
        experiments = flatten_experiments(params.input.experiments)
        if len(experiments) == 0:
            self.parser.print_help()
            return

        predicted_all = flex.reflection_table()

        for i_expt, expt in enumerate(experiments):
            if params.buffer_size > 0:
                # Hack to make the predicter predict reflections outside of the range
                # of the scan
                scan = expt.scan
                image_range = scan.get_image_range()
                oscillation = scan.get_oscillation()
                scan.set_image_range(
                    (
                        image_range[0] - params.buffer_size,
                        image_range[1] + params.buffer_size,
                    )
                )
                scan.set_oscillation(
                    (
                        oscillation[0] - params.buffer_size * oscillation[1],
                        oscillation[1],
                    )
                )

            # Populate the reflection table with predictions
            predicted = flex.reflection_table.from_predictions(
                expt, force_static=params.force_static, dmin=params.d_min
            )
            predicted["id"] = flex.int(len(predicted), i_expt)
            predicted_all.extend(predicted)

        # if we are not ignoring shadows, look for reflections in the masked
        # region, see https://github.com/dials/dials/issues/349

        if not params.ignore_shadows:
            shadowed = filter_shadowed_reflections(
                experiments, predicted_all, experiment_goniometer=True
            )
            predicted_all = predicted_all.select(~shadowed)

        try:
            predicted_all.compute_bbox(experiments)
        except Exception:
            pass

        # Save the reflections to file
        Command.start(f"Saving {len(predicted_all)} reflections to {params.output}")
        predicted_all.as_file(params.output)
        Command.end(f"Saved {len(predicted_all)} reflections to {params.output}")


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
