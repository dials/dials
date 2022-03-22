# LIBTBX_SET_DISPATCHER_NAME dials.filter_reflections


from __future__ import annotations

import logging
import token
from operator import itemgetter
from tokenize import TokenError, generate_tokens, untokenize

from cctbx import uctbx
from libtbx.phil import parse

from dials.algorithms.integration import filtering
from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors, tabulate
from dials.util.filter_reflections import SumAndPrfIntensityReducer, SumIntensityReducer
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

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

phil_scope = parse(phil_str, process_includes=True)


def eval_flag_expression(expression, reflections):
    """Test a Boolean expression of reflection flags for validity then
    evaluate it"""

    result = []
    tokens = iter(expression.split("\n"))
    g = generate_tokens(lambda: next(tokens))

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
            raise Sorry(f"errors found in {expression}")
        except StopIteration:
            break

        # Skip newline characters
        if toknum == token.NEWLINE:
            continue

        # Catch unwanted token types
        if toknum not in [token.OP, token.NAME, token.ENDMARKER]:
            raise Sorry(f"invalid token {token.tok_name[toknum]} found in {expression}")

        # Catch unwanted operators
        if toknum is token.OP and tokval not in "()|&~":
            raise Sorry(f"unrecognised operators found in {expression}")

        # Catch unrecognised flag names
        if toknum is token.NAME and tokval not in flag_names:
            raise Sorry(f"unrecognised flag name: {tokval}")

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
        print(tabulate(rows, header))
    else:
        print("No flags set")


def run_filtering(params, experiments, reflections):
    """Execute the script."""

    # Check params
    if params.d_min is not None and params.d_max is not None:
        if params.d_min > params.d_max:
            raise Sorry("d_min must be less than d_max")
    if params.d_min is not None or params.d_max is not None or params.ice_rings.filter:
        if "d" not in reflections:
            if experiments and any(experiments.crystals()):
                # Calculate d-spacings from the miller indices
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
            elif experiments:
                # Calculate d-spacings from the observed reflection centroids
                if "rlp" not in reflections:
                    if "xyzobs.mm.value" not in reflections:
                        reflections.centroid_px_to_mm(experiments)
                    reflections.map_centroids_to_reciprocal_space(experiments)
                d_star_sq = flex.pow2(reflections["rlp"].norms())
                reflections["d"] = uctbx.d_star_sq_as_d(d_star_sq)
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

    print(f"{len(reflections)} reflections loaded")

    # Filter by logical expression using flags
    if params.flag_expression is not None:
        inc = eval_flag_expression(params.flag_expression, reflections)
        reflections = reflections.select(inc)

    print(f"Selected {len(reflections)} reflections by flags")

    # Filter based on experiment ID
    if params.id:
        selection = reflections["id"] == params.id[0]
        for exp_id in params.id[1:]:
            selection = selection | (reflections["id"] == exp_id)
        reflections = reflections.select(selection)
        print(f"Selected {len(reflections)} reflections by experiment id")

    # Filter based on panel number
    if params.panel:
        selection = reflections["panel"] == params.panel[0]
        for pnl_id in params.panel[1:]:
            selection = selection | (reflections["panel"] == pnl_id)
        reflections = reflections.select(selection)
        print(f"Selected {len(reflections)} reflections by panel number")

    # Filter based on resolution
    if params.d_min is not None:
        selection = reflections["d"] >= params.d_min
        reflections = reflections.select(selection)
        print(f"Selected {len(reflections)} reflections with d >= {params.d_min:f}")

    # Filter based on resolution
    if params.d_max is not None:
        selection = reflections["d"] <= params.d_max
        reflections = reflections.select(selection)
        print(f"Selected {len(reflections)} reflections with d <= {params.d_max:f}")

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
        d_min = params.ice_rings.d_min
        width = params.ice_rings.width

        if d_min is None:
            d_min = flex.min(reflections["d"])

        ice_filter = filtering.PowderRingFilter(
            params.ice_rings.unit_cell,
            params.ice_rings.space_group.group(),
            d_min,
            width,
        )

        ice_sel = ice_filter(reflections["d"])

        print("Rejecting %i reflections at ice ring resolution" % ice_sel.count(True))
        reflections = reflections.select(~ice_sel)

    # Save filtered reflections to file
    if params.output.reflections:
        print(f"Saving {len(reflections)} reflections to {params.output.reflections}")
        reflections.as_file(params.output.reflections)


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


@show_mail_handle_errors()
def run(args=None):
    """Run the command line filtering script."""

    flags = list(flex.reflection_table.flags.names.items())
    flags.sort(key=itemgetter(0))

    # The script usage
    usage = "usage: dials.filter_reflections [options] experiment.expt"

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    log.config(verbosity=options.verbose)

    if not reflections:
        parser.print_help()
        return
    if len(reflections) != 1:
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
    run()
