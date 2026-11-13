# LIBTBX_SET_DISPATCHER_NAME dials.lambda_norm
"""
Normalise Laue/time-of-flight intensities as a function of wavelength.

A smooth w(λ) is fitted by Gaussian-process regression to leave-one-out ratios
of symmetry-equivalent intensities, and divided out of the profile and
summation intensities. Nothing else is corrected: run dials.scale afterwards.

The fitted curve is stored per reflection as lambda_norm.scale, with its
relative uncertainty as lambda_norm.sigma.

Examples::

  dials.lambda_norm integrated.expt integrated.refl

  dials.lambda_norm integrated.expt integrated.refl wavelength_range=2.1,3.9
"""

from __future__ import annotations

import json
import logging
import sys

from jinja2 import ChoiceLoader, Environment, PackageLoader

from dials.algorithms.scaling.tof.lambda_norm import (
    SCALE_COLUMN,
    WAVELENGTH,
    apply_lambda_norm,
    ascii_spectrum,
    fit_lambda_norm,
    lambda_norm_dict,
    lambda_norm_plots,
    normalisation_table,
    select_observations,
)
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.phil import parse
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.lambda_norm")

phil_scope = parse(
    """
output {
  experiments = lambda_norm.expt
    .type = str
  reflections = lambda_norm.refl
    .type = str
  log = dials.lambda_norm.log
    .type = str
  json = None
    .type = str
    .help = "Write the fitted curve and its training data to this json file."
  html = None
    .type = str
    .help = "Write an html report of the fit to this file."
}
intensity_choice = *auto profile sum
  .type = choice
  .help = "The intensities the fit is trained on: auto takes profile-fitted"
          "intensities if present. Profile and summation intensities are both"
          "corrected."
min_isigi = -5.0
  .type = float
  .help = "Observations below this I/σ are not used in the fit. Weak and negative"
          "intensities are kept by default: cutting them biases w(λ) at the"
          "ends of the spectrum."
partiality_cutoff = 0.25
  .type = float
  .help = "Observations below this partiality are not used in the fit, and"
          "summation intensities are divided by partiality, as dials.scale does."
wavelength_range = None
  .type = floats(size=2, value_min=0)
  .help = "Fit over this wavelength range (Å); w(λ) is held at its end values"
          "outside it. Default: the range of the data."
reference_wavelength = None
  .type = float(value_min=0)
  .help = "Wavelength (Å) at which w(λ) = 1. Default: the median wavelength"
          "of the observations."
anomalous = False
  .type = bool
  .help = "Keep Friedel mates separate when grouping equivalents."
gpr {
  nbins = 50
    .type = int(value_min=4)
    .help = "Number of equal-count wavelength bins the GP is trained on. Reduce"
            "it if the length scale reaches its lower search limit."
  lengthscale = None
    .type = float(value_min=0)
    .help = "Fix the GP length scale (Å at the reference wavelength) instead of"
            "choosing it by marginal likelihood."
  kernel = *squared_exponential matern32
    .type = choice
    .help = "Matérn-3/2 is less smooth; use it if the true curve has kinks."
  n_cycles = 3
    .type = int(value_min=1)
    .help = "Number of fit cycles. Each cycle divides the current w(λ) out of"
            "the equivalents before forming the ratios, so the ratios do not"
            "depend on the wavelengths the equivalents happened to be measured at."
}
weighting = *expected measured
  .type = choice
  .help = "expected: weight each ratio by the variance the observation would"
          "have at its predicted intensity. With counting statistics the"
          "measured σ is small exactly when I came out low, which turns the"
          "bin average into a harmonic mean."
per_experiment_factor = True
  .type = bool
  .help = "Divide a constant factor per experiment out of the ratios during the"
          "fit, standing in for the per-run scale that dials.scale refines"
          "afterwards. The factors are not applied to the output."
"""
)


@show_mail_handle_errors()
def run(args: list[str] | None = None, phil=phil_scope) -> None:
    usage = "dials.lambda_norm integrated.expt integrated.refl [options]"
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=False)
    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if len(reflections) != 1:
        raise Sorry("dials.lambda_norm takes a single reflection file")
    reflections = reflections[0]
    if WAVELENGTH not in reflections:
        raise Sorry(f"No {WAVELENGTH} column: is this time-of-flight data?")
    if SCALE_COLUMN in reflections:
        raise Sorry("These reflections are already wavelength-normalised")

    try:
        observations = select_observations(reflections, experiments, params)
        logger.info(
            "Fitting w(λ) to %d observations in %d groups of equivalents",
            len(observations.intensity),
            len(observations.group_offsets) - 1,
        )
        result = fit_lambda_norm(observations, params)
    except ValueError as e:
        raise Sorry(e)

    logger.info(
        "\nWavelength normalisation, w(λ) = 1 at %.3f Å:\n%s\n\n%s",
        result.gpr.lambda_ref(),
        normalisation_table(result.gpr),
        ascii_spectrum(result.gpr),
    )
    reflections = apply_lambda_norm(reflections, result.gpr)

    logger.info("Saving reflections to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)
    logger.info("Saving experiments to %s", params.output.experiments)
    experiments.as_file(params.output.experiments)

    if params.output.json or params.output.html:
        data = lambda_norm_dict(result)
    if params.output.json:
        logger.info("Writing json data to %s", params.output.json)
        with open(params.output.json, "w") as f:
            json.dump(data, f, indent=2)
    if params.output.html:
        logger.info("Writing html report to %s", params.output.html)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        html = (
            Environment(loader=loader)
            .get_template("simple_report.html")
            .render(
                page_title="dials.lambda_norm report",
                panel_title="Wavelength normalisation",
                panel_id="lambda_norm",
                graphs=lambda_norm_plots(data),
            )
        )
        with open(params.output.html, "wb") as f:
            f.write(html.encode("utf-8", "xmlcharrefreplace"))


if __name__ == "__main__":
    run()
