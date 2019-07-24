#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dials.filter_reflections

from __future__ import absolute_import, division, print_function
import logging
from operator import itemgetter
import token
from tokenize import generate_tokens, TokenError, untokenize
from StringIO import StringIO

from dials.util import Sorry, log, show_mail_on_error
from dials.util.filter_reflections import SumAndPrfIntensityReducer, SumIntensityReducer
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dials.algorithms.integration import filtering
from dials.algorithms.spot_finding.per_image_analysis import map_to_reciprocal_space
from libtbx.phil import parse
from libtbx.table_utils import simple_table
from cctbx import uctbx


logger = logging.getLogger("dials")

help_message = """

This program takes reflection files as input and filters them based on user-
specified criteria, to write out a subset of the original file.

Filtering is first done by evaluating the optional boolean 'flag_expression'
using reflection flag values. The operators allowed are '&' for 'and', '|' for
'or', and '~' for 'not'. Expressions may contain nested sub-expressions using
parentheses.

Following this, optional additional filters are applied according to values in
the reflection table, such as by resolution or user-defined masks.

If a reflection file is passed in to the program but no filtering parameters
are set, a table will be printed, giving the flag values present in the
reflection file.


Examples::

  dials.filter_reflections refined.refl \
    flag_expression=used_in_refinement

  dials.filter_reflections integrated.refl \
    flag_expression="integrated & ~reference_spot"

  dials.filter_reflections integrated.refl \
    flag_expression="indexed & (failed_during_summation | failed_during_profile_fitting)"

  dials.filter_reflections indexed.refl indexed.expt \
    d_max=20 d_min=2.5
"""

phil_str = """
    output {
    reflections = 'filtered.refl'
        .type = str
        .help = "The filtered reflections output filename"
    }

    flag_expression = None
    .type = str
    .help = "Boolean expression to select reflections based on flag values"

    id = None
    .type = ints(value_min=0)
    .help = "Select reflections by experiment IDs"

    panel = None
    .type = ints(value_min=0)
    .help = "Select reflections by panels they intersect"

    d_min = None
    .type = float
    .help = "The maximum resolution"

    d_max = None
    .type = float
    .help = "The minimum resolution"

    partiality {
    min = None
        .type = float(value_min = 0, value_max = 1)
        .help = "The minimum reflection partiality for inclusion."
    max = None
        .type = float(value_min = 0, value_max = 1)
        .help = "The maximum reflection partiality for inclusion."
    }

    select_good_intensities = False
    .type = bool
    .help = "Combined filter to select only fully integrated and"
            "trustworthy intensities"

    dead_time
    {
    value = 0
        .help = "Detector dead time in ms, assumed to be at the end of the"
                "exposure time."
        .type = float(value_min=0)

    reject_fraction = 0
        .help = "Reject reflections which overlap by more than the given"
                "fraction with the dead region of the image."
        .type = float(value_min=0, value_max=1)
    }

    include scope dials.util.masking.phil_scope ice_rings

"""


def eval_flag_expression(expression, reflections):
    """Test a Boolean expression of reflection flags for validity then
    evaluate it"""

    result = []
    g = generate_tokens(StringIO(expression).readline)

    flags = list(flex.reflection_table.flags.names.items())
    flags.sort(key=itemgetter(0))
    flag_names, _ = zip(*flags)

    # define shorthand function
    def get_flag(flag):
        return reflections.get_flags(getattr(reflections.flags, flag))

    while True:

        # Extract next token, catching unmatched brackets
        try:
            toknum, tokval, _, _, _ = next(g)
        except TokenError:
            raise Sorry("errors found in {}".format(expression))
        except StopIteration:
            break

        # Catch unwanted token types
        if toknum not in [token.OP, token.NAME, token.ENDMARKER]:
            raise Sorry("invalid tokens found in {}".format(expression))

        # Catch unwanted operators
        if toknum is token.OP and tokval not in "()|&~":
            raise Sorry("unrecognised operators found in {}".format(expression))

        # Catch unrecognised flag names
        if toknum is token.NAME and tokval not in flag_names:
            raise Sorry("unrecognised flag name: {}".format(tokval))

        # Replace names with valid lookups in the reflection table
        if toknum is token.NAME:
            # ("NAME", "get_flags")
            # ("OP", "(")
            # ("STRING", "'indexed'")
            # ("OP", ")")
            result.extend(
                [
                    (token.NAME, "get_flag"),
                    (token.OP, "("),
                    (token.STRING, repr(tokval)),
                    (token.OP, ")"),
                ]
            )
        else:
            result.append((toknum, tokval))

    # Evaluate and return the result
    return eval(untokenize(result), {"get_flag": get_flag, "__builtins__": None}, {})


def run_analysis(flags, reflections):
    """Print a table of flags present in the reflections file"""

    header = ["flag", "nref"]
    rows = []
    for name, val in flags:
        n = (reflections.get_flags(val)).count(True)
        if n > 0:
            rows.append([name, "%d" % n])
    if rows:
        st = simple_table(rows, header)
        print(st.format())
    else:
        print("No flags set")

    return


def run_filtering(params, experiments, reflections):
    """Execute the script."""

    # Check params
    if params.d_min is not None and params.d_max is not None:
        if params.d_min > params.d_max:
            raise Sorry("d_min must be less than d_max")
    if params.d_min is not None or params.d_max is not None:
        if "d" not in reflections:
            if experiments:
                print(
                    "Reflection table does not have resolution information. "
                    "Attempting to calculate this from the experiment list"
                )
                sel = reflections["id"] >= 0
                if sel.count(False) > 0:
                    print(
                        "Removing {} reflections with negative experiment id".format(
                            sel.count(False)
                        )
                    )
                reflections = reflections.select(sel)
                reflections.compute_d(experiments)
            else:
                raise Sorry(
                    "reflection table has no resolution information "
                    "and no experiment list provided to calculate it"
                )

    # Check params
    if params.partiality.min is not None and params.partiality.max is not None:
        if params.min > params.max:
            raise Sorry("partiality.min must be less than partiality.d_max")
    if params.partiality.min is not None or params.partiality.max is not None:
        if "partiality" not in reflections:
            raise Sorry("Reflection table has no partiality information")

    print("{} reflections loaded".format(len(reflections)))

    # Filter by logical expression using flags
    if params.flag_expression is not None:
        inc = eval_flag_expression(params.flag_expression, reflections)
        reflections = reflections.select(inc)

    print("Selected {} reflections by flags".format(len(reflections)))

    # Filter based on experiment ID
    if params.id:
        selection = reflections["id"] == params.id[0]
        for exp_id in params.id[1:]:
            selection = selection | (reflections["id"] == exp_id)
        reflections = reflections.select(selection)
        print("Selected %d reflections by experiment id" % (len(reflections)))

    # Filter based on panel number
    if params.panel:
        selection = reflections["panel"] == params.panel[0]
        for pnl_id in params.panel[1:]:
            selection = selection | (reflections["panel"] == pnl_id)
        reflections = reflections.select(selection)
        print("Selected %d reflections by panel number" % (len(reflections)))

    # Filter based on resolution
    if params.d_min is not None:
        selection = reflections["d"] >= params.d_min
        reflections = reflections.select(selection)
        print("Selected %d reflections with d >= %f" % (len(reflections), params.d_min))

    # Filter based on resolution
    if params.d_max is not None:
        selection = reflections["d"] <= params.d_max
        reflections = reflections.select(selection)
        print("Selected %d reflections with d <= %f" % (len(reflections), params.d_max))

    # Filter based on partiality
    if params.partiality.min is not None:
        selection = reflections["partiality"] >= params.partiality.min
        reflections = reflections.select(selection)
        print(
            "Selected %d reflections with partiality >= %f"
            % (len(reflections), params.partiality.min)
        )

    # Filter based on partiality
    if params.partiality.max is not None:
        selection = reflections["partiality"] <= params.partiality.max
        reflections = reflections.select(selection)
        print(
            "Selected %d reflections with partiality <= %f"
            % (len(reflections), params.partiality.max)
        )

    # 'Good' intensity selection
    if params.select_good_intensities:
        if "intensity.sum.variance" not in reflections:
            raise Sorry("No intensity.sum.variance in reflection table.")
        if "intensity.prf.variance" in reflections:
            reducer = SumAndPrfIntensityReducer
        else:
            reducer = SumIntensityReducer
        try:
            reflections = reducer.filter_bad_variances(reflections)
            reflections = reducer.combine_and_filter_partials(
                reflections, partiality_threshold=0.99, combine_partials=False
            )
        except ValueError as e:
            raise Sorry(e)

    # Dead time filter
    if params.dead_time.value > 0:
        reflections = filter_by_dead_time(
            reflections,
            experiments,
            params.dead_time.value,
            params.dead_time.reject_fraction,
        )

    # Filter powder rings
    if params.ice_rings.filter:
        if "d" in reflections:
            d_spacings = reflections["d"]
        else:
            if "rlp" not in reflections:
                imageset = experiments.imagesets()[0]
                assert imageset is not None
                reflections = map_to_reciprocal_space(reflections, imageset)
            d_star_sq = flex.pow2(reflections["rlp"].norms())
            d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

        d_min = params.ice_rings.d_min
        width = params.ice_rings.width

        if d_min is None:
            d_min = flex.min(d_spacings)

        ice_filter = filtering.PowderRingFilter(
            params.ice_rings.unit_cell,
            params.ice_rings.space_group.group(),
            d_min,
            width,
        )

        ice_sel = ice_filter(d_spacings)

        print("Rejecting %i reflections at ice ring resolution" % ice_sel.count(True))
        reflections = reflections.select(~ice_sel)

    # Save filtered reflections to file
    if params.output.reflections:
        print(
            "Saving {0} reflections to {1}".format(
                len(reflections), params.output.reflections
            )
        )
        reflections.as_pickle(params.output.reflections)

    return


def filter_by_dead_time(reflections, experiments, dead_time=0, reject_fraction=0):
    if not experiments:
        raise Sorry("an experiment list must be provided to filter by dead time")

    if len(experiments) > 1:
        raise Sorry("only a single experiment is supported in this mode")
    experiment = experiments[0]

    sel = reflections.get_flags(reflections.flags.integrated)
    reflections = reflections.select(sel)

    if not reflections:
        raise Sorry("no integrated reflections present")

    goniometer = experiment.goniometer
    beam = experiment.beam

    phi_range = reflections.compute_phi_range(
        goniometer.get_rotation_axis(),
        beam.get_s0(),
        experiment.profile.sigma_m(deg=False),
        experiment.profile.n_sigma(),
    )
    phi1, phi2 = phi_range.parts()

    scan = experiment.scan
    exposure_time = scan.get_exposure_times()[0]
    assert scan.get_exposure_times().all_eq(exposure_time)
    phi_start, phi_width = scan.get_oscillation(deg=False)
    phi_range_dead = phi_width * (dead_time / 1000) / exposure_time

    sel_good = flex.bool(len(reflections), True)

    start, end = scan.get_array_range()
    for i in range(start, end):
        phi_dead_start = phi_start + (i + 1) * phi_width - phi_range_dead
        phi_dead_end = phi_dead_start + phi_range_dead

        left = phi1.deep_copy()
        left.set_selected(left < phi_dead_start, phi_dead_start)

        right = phi2.deep_copy()
        right.set_selected(right > phi_dead_end, phi_dead_end)

        overlap = (right - left) / (phi2 - phi1)

        sel = overlap > reject_fraction

        sel_good.set_selected(sel, False)
        print("Rejecting %i reflections from image %i" % (sel.count(True), i))

    print(
        "Keeping %i reflections (rejected %i)"
        % (sel_good.count(True), sel_good.count(False))
    )

    return reflections.select(sel_good)


def run():
    """Run the command line filtering script."""

    flags = list(flex.reflection_table.flags.names.items())
    flags.sort(key=itemgetter(0))

    phil_scope = parse(phil_str, process_includes=True)

    # The script usage
    usage = "usage: dials.filter_reflections [options] experiment.expt"

    # Create the parser
    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    log.config(verbosity=options.verbose)

    if not reflections:
        parser.print_help()
        raise Sorry("No valid reflection file given")
    if len(reflections) != 1:
        parser.print_help()
        raise Sorry("Exactly 1 reflection file must be specified")
    reflections = reflections[0]

    # Check if any filter has been set using diff_phil
    filter_def = [
        o for o in parser.diff_phil.objects if o.name not in ["input", "output"]
    ]
    if not filter_def:
        print("No filter specified. Performing analysis instead.")
        run_analysis(flags, reflections)
    else:
        run_filtering(params, experiments, reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        run()
