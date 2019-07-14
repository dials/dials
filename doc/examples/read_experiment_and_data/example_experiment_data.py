#!/usr/bin/env dials.python
#
# example_experiment_data.py
#
#  Copyright (C) 2017 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# Example code for how to load experiments and reflections in the DIALS
# framework

from __future__ import absolute_import, division, print_function
from libtbx.phil import parse

help_message = """

pass experiment.expt indexed.refl

"""

phil_scope = parse(
    """
png = 'example.png'
  .type = str
  .help = 'Output name for .png'
""",
    process_includes=True,
)


class Script(object):
    """A class for running the script."""

    def __init__(self):
        from dials.util.options import OptionParser
        import libtbx.load_env

        usage = (
            "usage: %s [options] indexed.expt indexed.refl" % libtbx.env.dispatcher_name
        )

        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            check_format=False,
            read_experiments=True,
            read_reflections=True,
        )

    def run(self):
        from dials.array_family import flex  # import dependency
        from scitbx import matrix
        from dials.util.options import flatten_experiments
        from dials.util.options import flatten_reflections

        params, options = self.parser.parse_args(show_diff_phil=True)

        experiments = flatten_experiments(params.input.experiments)
        reflections = flatten_reflections(params.input.reflections)

        if len(experiments) == 0:
            self.parser.print_help()
            return

        if len(reflections) != 1:
            self.parser.print_help()
            return

        reflections = reflections[0]

        print("Read %d reflections" % len(reflections))

        indexed = reflections.select(reflections.get_flags(reflections.flags.indexed))

        print("Kept %d indexed reflections" % len(indexed))

        for name in sorted(indexed.keys()):
            print("Found column %s" % name)

        for reflection in indexed[:3]:
            print(reflection)

        # verify that these experiments correspond to exactly one imageset, one
        # detector, one beam (obviously)
        for experiment in experiments[1:]:
            assert experiment.imageset == experiments[0].imageset
            assert experiment.beam == experiments[0].beam
            assert experiment.detector == experiments[0].detector

        # now perform some calculations - the only things different from one
        # experiment to the next will be crystal models
        crystals = [experiment.crystal for experiment in experiments]
        detector = experiments[0].detector
        beam = experiments[0].beam
        imageset = experiments[0].imageset

        # derived quantities
        wavelength = beam.get_wavelength()
        s0 = matrix.col(beam.get_s0())

        # in here do some jiggery-pokery to cope with this being interpreted as
        # a rotation image in here i.e. if scan is not None; derive goniometer
        # matrix else set this to identity

        scan = experiments[0].scan
        goniometer = experiments[0].goniometer

        if scan and goniometer:
            angle = scan.get_angle_from_array_index(
                0.5 * sum(imageset.get_array_range())
            )
            axis = matrix.col(goniometer.get_rotation_axis_datum())
            F = matrix.sqr(goniometer.get_fixed_rotation())
            S = matrix.sqr(goniometer.get_setting_rotation())
            R = S * axis.axis_and_angle_as_r3_rotation_matrix(angle, deg=True) * F
        else:
            R = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

        assert len(detector) == 1


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
