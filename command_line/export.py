from __future__ import absolute_import, division, print_function

import logging
import sys

from libtbx.phil import parse
from libtbx import Auto

logger = logging.getLogger("dials.command_line.export")

help_message = """

This program is used to export the results of dials processing in various
formats.

The output formats currently supported are:

MTZ format exports the files as an unmerged mtz file, ready for input to
downstream programs such as Pointless and Aimless. For exporting integrated,
but unscaled data, the required input is a models.expt file and an
integrated.refl file. For exporting scaled data, the required input is a
models.expt file and a scaled.pickle file, also passing the option
intensity=scale.

NXS format exports the files as an NXmx file. The required input is a
models.expt file and an integrated.refl file.

MMCIF format exports the files as an mmcif file. The required input is a
models.expt file and an integrated.refl file.

XDS_ASCII format exports intensity data and the experiment metadata in the
same format as used by the output of XDS in the CORRECT step - output can
be scaled with XSCALE.

SADABS format exports intensity data (and geometry by direction cosines)
as an ersatz-SADABS format reverse engineered from the file format used by
EvalCCD for input to SADABS.

MOSFLM format exports the files as an index.mat mosflm-format matrix file and a
mosflm.in file containing basic instructions for input to mosflm. The required
input is an models.expt file.

XDS format exports a models.expt file as XDS.INP and XPARM.XDS files. If a
reflection pickle is given it will be exported as a SPOT.XDS file.

Examples::

  # Export to mtz
  dials.export models.expt integrated.refl
  dials.export models.expt integrated.refl mtz.hklout=integrated.mtz
  dials.export models.expt scaled.pickle intensity=scale mtz.hklout=scaled.mtz

  # Export to nexus
  dials.export models.expt integrated.refl format=nxs
  dials.export models.expt integrated.refl format=nxs nxs.hklout=integrated.nxs

  # Export to mmcif
  dials.export models.expt integrated.refl format=mmcif
  dials.export models.expt integrated.refl format=mmcif mmcif.hklout=integrated.mmcif

  # Export to mosflm
  dials.export models.expt integrated.refl format=mosflm

  # Export to xds
  dials.export strong.pickle format=xds
  dials.export indexed.pickle format=xds
  dials.export models.expt format=xds
  dials.export models.expt indexed.pickle format=xds
"""

phil_scope = parse(
    """

  format = *mtz sadabs nxs mmcif mosflm xds xds_ascii json
    .type = choice
    .help = "The output file format"

  intensity = *auto profile sum scale
    .type = choice(multi=True)
    .help = "Choice of which intensities to export. Allowed combinations:
            scale, profile, sum, profile+sum, sum+profile+scale. Auto will
            default to scale or profile+sum depending on if the data are scaled."


  debug = False
    .type = bool
    .help = "Output additional debugging information"

  mtz {

    combine_partials = True
      .type = bool
      .help = "Combine partials that have the same partial id into one
        reflection, with an updated partiality given by the sum of the
        individual partialities."

    partiality_threshold=0.99
      .type = float
      .help = "All reflections with partiality values above the partiality
        threshold will be retained. This is done after any combination of
        partials if applicable."

    min_isigi = -5
      .type = float
      .help = "Exclude reflections with unfeasible values of I/Sig(I)"

    force_static_model = False
      .type = bool
      .help = "Force program to use static model even if scan varying is present"

    filter_ice_rings = False
      .type = bool
      .help = "Filter reflections at ice ring resolutions"

    d_min = None
      .type = float
      .help = "Filter out reflections with d-spacing below d_min"

    hklout = auto
      .type = path
      .help = "The output MTZ filename, defaults to integrated.mtz or scaled_unmerged.mtz"
              "depending on if the input data are scaled."

    crystal_name = XTAL
      .type = str
      .help = "The name of the crystal, for the mtz file metadata"
  }

  sadabs {

    hklout = integrated.sad
      .type = path
      .help = "The output raw sadabs file"
    run = 1
      .type = int
      .help = "Batch number / run number for output file"
    predict = False
      .type = bool
      .help = "Compute centroids with static model, not observations"

  }

  xds_ascii {

    hklout = DIALS.HKL
      .type = path
      .help = "The output raw hkl file"

  }

  nxs {

    hklout = integrated.nxs
      .type = path
      .help = "The output Nexus file"

  }

  mmcif {

    hklout = auto
      .type = path
      .help = "The output CIF file, defaults to integrated.cif or scaled_unmerged.cif
        depending on if the data are scaled."

  }

  mosflm {

    directory = mosflm
      .type = path
      .help = "The output directory for mosflm output"

  }

  xds {

    directory = xds
      .type = path
      .help = "The output directory for xds output"

  }

  json {
    filename = rlp.json
      .type = path
    compact = True
      .type = bool
    n_digits = 6
      .type = int(value_min=1)
      .help = "Number of decimal places to be used for representing the"
              "reciprocal lattice points."
  }

  output {
    log = dials.export.log
      .type = path
      .help = "The log filename"
  }
"""
)


def _check_input(experiments, reflections, params=None, check_intensities=False):
    if not experiments:
        sys.exit("export requires an experiment list")
    if len(reflections) != 1:
        sys.exit("export requires 1 reflection table")
    if check_intensities:
        if not params:
            sys.exit("No parameters given")
        if "profile" not in params.intensity and "sum" not in params.intensity:
            sys.exit(
                "Only intensity options containing sum or profile can be exported in this format"
            )
        if (
            "intensity.sum.value" not in reflections[0]
            and "intensity.prf.value" not in reflections[0]
        ):
            sys.exit(
                "Unable to find 'intensity.sum.value' or 'intensity.prf.value' columns in reflection table."
            )


def export_mtz(params, experiments, reflections):
    """
    Export reflections in MTZ format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, reflections)

    from dials.util.export_mtz import export_mtz

    # Handle case where user has passed data before integration
    if (
        "intensity.sum.value" not in reflections[0]
        or "intensity.prf.value" not in reflections[0]
    ):
        sys.exit(
            "Error: No intensity data in reflections; cannot export un-integrated data to MTZ"
        )

    try:
        m = export_mtz(reflections[0], experiments, params)
    except ValueError as e:
        sys.exit(e)
    from six.moves import cStringIO as StringIO

    summary = StringIO()
    m.show_summary(out=summary)
    logger.info("")
    logger.info(summary.getvalue())


def export_sadabs(params, experiments, reflections):
    """
    Export data in HKL format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, reflections, params=params, check_intensities=True)

    from dials.util.export_sadabs import export_sadabs

    try:
        export_sadabs(reflections[0], experiments, params)
    except ValueError as e:
        sys.exit(e)


def export_xdsascii(params, experiments, reflections):
    """
    Export data in XDS_ASCII format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, reflections, params=params, check_intensities=True)

    from dials.util.export_xds_ascii import export_xds_ascii

    try:
        export_xds_ascii(reflections[0], experiments, params)
    except ValueError as e:
        sys.exit(e)


def export_nexus(params, experiments, reflections):
    """
    Export data in Nexus format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, reflections)

    from dials.util.nexus import dump

    dump(experiments, reflections[0], params.nxs.hklout)


def export_mmcif(params, experiments, reflections):
    """
    Export data in CIF format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, reflections)

    from dials.util.export_mmcif import MMCIFOutputFile

    outfile = MMCIFOutputFile(params)
    try:
        outfile.write(experiments, reflections[0])
    except ValueError as e:
        sys.exit(e)


def export_mosflm(params, experiments, reflections):
    """
    Export stuff in mosflm format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    _check_input(experiments, [None])

    if reflections:
        sys.exit("Mosflm exporter does not need a reflection table")

    from dials.util.mosflm import dump

    dump(experiments, params.mosflm.directory)


def export_xds(params, experiments, reflections):
    """
    Export stuff in xds format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    # Check the input
    if len(reflections) > 1:
        sys.exit("XDS exporter requires 0 or 1 reflection table")

    if reflections:
        reflections = reflections[0]

    from dials.util.xds import dump

    dump(experiments, reflections, params.xds.directory)


def export_json(params, experiments, reflections):
    """
    Export reflections in json format

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables
    """

    # Check the input
    _check_input(experiments, [None])
    if not reflections:
        sys.exit("json exporter requires a reflection table")

    from dials.util import export_json
    from scitbx.array_family import flex

    imagesets = [expt.imageset for expt in experiments]

    if len(reflections) != len(imagesets):
        sys.exit(
            "Mismatch between %d reflections lists and %d imagesets"
            % (len(reflections), len(imagesets))
        )
    selected_reflections = None
    for i, refl in enumerate(reflections):
        refl["imageset_id"] = flex.size_t(refl.size(), i)
        if selected_reflections is None:
            selected_reflections = refl
        else:
            selected_reflections.extend(refl)

    exporter = export_json.ReciprocalLatticeJson(experiments, selected_reflections)
    exporter.as_json(
        filename=params.json.filename,
        compact=params.json.compact,
        n_digits=params.json.n_digits,
        experiments=experiments,
    )


if __name__ == "__main__":
    from dials.util.options import OptionParser, reflections_and_experiments_from_files
    from dials.util.version import dials_version
    from dials.util import log

    usage = "dials.export models.expt reflections.pickle [options]"

    # Create the option parser
    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        phil=phil_scope,
        epilog=help_message,
    )

    # Get the parameters
    params, options = parser.parse_args(show_diff_phil=False)

    # Configure the logging
    log.config(logfile=params.output.log)

    # Print the version number
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if not params.input.experiments and not params.input.reflections:
        parser.print_help()
        sys.exit()

    # Get the experiments and reflections
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # do auto intepreting of intensity choice:
    # note that this may still fail certain checks further down the processing,
    # but these are the defaults to try
    if params.intensity in ([None], [Auto], ["auto"]) and reflections:
        if ("intensity.scale.value" in reflections[0]) and (
            "intensity.scale.variance" in reflections[0]
        ):
            params.intensity = ["scale"]
            logger.info("Data appears to be scaled, setting intensity = scale")
        else:
            params.intensity = ["profile", "sum"]
            logger.info("Data appears to be unscaled, setting intensity = profile+sum")

    # Choose the exporter
    exporter = {
        "mtz": export_mtz,
        "sadabs": export_sadabs,
        "xds_ascii": export_xdsascii,
        "nxs": export_nexus,
        "mmcif": export_mmcif,
        "mosflm": export_mosflm,
        "xds": export_xds,
        "json": export_json,
    }.get(params.format)
    if not exporter:
        sys.exit("Unknown format: %s" % params.format)

    # Export the data
    exporter(params, experiments, reflections)
