"""
This program imports xds processed data for use in dials.

It requires up to three things to create and experiment list and reflection table.
    - an XDS.INP, to specify the geometry,
    - one of "XDS_ASCII.HKL", "INTEGRATE.HKL", "GXPARM.XDS", "XPARM.XDS", which is needed to create the experiment
    - INTEGRATE.HKL or SPOT.XDS file to create a reflection table.

To run the program, the easiest thing to do is provide a directory containing these files

Example use cases:

dials.import_xds /path/to/folder/containing/xds/inp/

dials.import_xds /path/to/folder/containing/xds/inp/INTEGRATE.HKL

dials.import_xds /path/to/folder/containing/xds/inp/ SPOT.XDS   # be explicit about which file to use to create reflections (default is to use INTEGRATE.HKL)

dials.import_xds /path/to/folder/containing/xds/inp/ xds_file=XPARM.XDS   # specify which extra file should be used to create experiment metadata

dials.import_xds /path/to/folder/containing/xds/inp/ /path/to/INTEGRATE.HKL   # will take XDS.INP from the directory, and everything else needed from the specified INTEGRATE.HKL file
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Optional

from cctbx import sgtbx
from dxtbx.model import Crystal
from dxtbx.model.experiment_list import ExperimentListFactory
from iotbx.xds import integrate_hkl, spot_xds
from libtbx.phil import parse
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix

from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.import_xds")

required_files_to_make_experiments = [
    "XDS_ASCII.HKL",
    "INTEGRATE.HKL",
    "GXPARM.XDS",
    "XPARM.XDS",
]


class SpotXDSImporter:
    """Class to import a spot.xds file to a reflection table."""

    def __init__(self, spot_xds):
        self._spot_xds = spot_xds

    def __call__(self, params, options):
        """Import the spot.xds file."""
        # Read the SPOT.XDS file
        logger.info(f"Reading {self._spot_xds}")
        handle = spot_xds.reader()
        handle.read_file(self._spot_xds)
        centroid = handle.centroid
        intensity = handle.intensity
        try:
            miller_index = handle.miller_index
        except AttributeError:
            miller_index = None
        logger.info(f"Read {len(centroid)} spots from {self._spot_xds}")

        # Create the reflection list
        logger.info("Creating reflection list")
        table = flex.reflection_table()
        table["id"] = flex.int(len(centroid), 0)
        table["panel"] = flex.size_t(len(centroid), 0)
        if miller_index:
            table["miller_index"] = flex.miller_index(miller_index)
        table["xyzobs.px.value"] = flex.vec3_double(centroid)
        table["intensity.sum.value"] = flex.double(intensity)

        # Remove invalid reflections
        logger.info("Removing invalid reflections")
        if miller_index and params.remove_invalid:
            flags = flex.bool([h != (0, 0, 0) for h in table["miller_index"]])
            table = table.select(flags)
        logger.info(f"Removed invalid reflections, {len(table)} remaining")

        # Fill empty standard columns
        if params.add_standard_columns:
            logger.info("Adding standard columns")
            rt = flex.reflection_table.empty_standard(len(table))
            rt.update(table)
            table = rt
            # set variances to unity
            table["xyzobs.mm.variance"] = flex.vec3_double(len(table), (1, 1, 1))
            table["xyzobs.px.variance"] = flex.vec3_double(len(table), (1, 1, 1))
            logger.info("Standard columns added")

        # Output the table to pickle file
        if params.output.reflections is None:
            params.output.reflections = "spot_xds.refl"
        logger.info(f"Saving reflection table to {params.output.reflections}")
        table.as_file(params.output.reflections)
        logger.info(f"Saved reflection table to {params.output.reflections}")


class IntegrateHKLImporter:
    """Class to import an INTEGRATE.HKL file to a reflection table."""

    def __init__(self, integrate_hkl, experiments):
        self._integrate_hkl = integrate_hkl
        self._experiment = experiments[0]
        self._experiments = experiments

    def __call__(self, params, options):
        """Import the integrate.hkl file."""
        # Get the unit cell to calculate the resolution
        uc = self._experiment.crystal.get_unit_cell()

        # Read the INTEGRATE.HKL file
        logger.info(f"Reading {self._integrate_hkl}")
        handle = integrate_hkl.reader()
        handle.read_file(self._integrate_hkl)
        hkl = flex.miller_index(handle.hkl)
        xyzcal = flex.vec3_double(handle.xyzcal)
        xyzobs = flex.vec3_double(handle.xyzobs)
        iobs = flex.double(handle.iobs)
        sigma = flex.double(handle.sigma)

        # for INTEGRATE.HKL this contains /only/ L - not P!
        rlp = flex.double(handle.rlp)

        peak = flex.double(handle.peak) * 0.01
        corr = flex.double(handle.corr) * 0.01
        if len(handle.iseg):
            panel = flex.size_t(handle.iseg) - 1
        else:
            panel = flex.size_t(len(hkl), 0)
        logger.info(f"Read {len(hkl)} reflections from {self._integrate_hkl}")

        undo_ab = True

        # ideally we would like to do that - but because of limited precision in
        # INTEGRATE.HKL we could get into trouble here ...
        if handle.variance_model and undo_ab:
            variance_model_a = handle.variance_model[0]
            variance_model_b = handle.variance_model[1]
            logger.info(
                "Undoing input variance model a, b = %s %s"
                % (variance_model_a, variance_model_b)
            )

            # undo variance model:
            vari0 = flex.pow2(sigma * peak / rlp)
            isq = flex.pow2(iobs * peak / rlp)
            vari1 = (vari0 / variance_model_a) - (variance_model_b * isq)
            for i in range(0, (sigma.size() - 1)):
                if vari1[i] <= 0:
                    logger.info("WARNING:", i, iobs[i], sigma[i], vari1[i])
                    vari1[i] = sigma[i] * sigma[i]
            sigma = flex.sqrt(vari1) * rlp / peak

        if len(self._experiment.detector) > 1:
            for p_id, p in enumerate(self._experiment.detector):
                sel = panel == p_id
                offset = p.get_raw_image_offset()
                xyzcal.set_selected(sel, xyzcal.select(sel) - (offset[0], offset[1], 0))
                xyzobs.set_selected(sel, xyzobs.select(sel) - (offset[0], offset[1], 0))

        # Derive the reindex matrix
        rdx = self.derive_reindex_matrix(handle)
        logger.info("Reindex matrix:\n%d %d %d\n%d %d %d\n%d %d %d" % (rdx.elems))

        # Reindex the reflections
        logger.info("Reindexing reflections")
        cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(sgtbx.rot_mx(rdx.elems)))
        hkl = cb_op.apply(hkl)
        logger.info(f"Reindexed {len(hkl)} reflections")

        # Create the reflection list
        logger.info("Creating reflection table")
        table = flex.reflection_table()
        table["id"] = flex.int(len(hkl), 0)
        table["panel"] = panel
        table["miller_index"] = hkl

        table["xyzcal.px"] = xyzcal
        table["xyzobs.px.value"] = xyzobs
        # uncorrected "partials":
        table["intensity.prf.value"] = iobs * peak / rlp
        table["intensity.prf.variance"] = flex.pow2(sigma * peak / rlp)
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)

        table["d"] = flex.double(uc.d(h) for h in hkl)

        table["partiality"] = peak
        table["profile.correlation"] = corr

        table.centroid_px_to_mm(self._experiments)

        # same as prf?
        table["intensity.sum.value"] = iobs * peak / rlp
        table["intensity.sum.variance"] = flex.pow2(sigma * peak / rlp)
        table.map_centroids_to_reciprocal_space(self._experiments)

        # LP-corrected fulls:
        table["intensity.cor.value"] = iobs
        table["intensity.cor.variance"] = flex.pow2(sigma)

        table["batch"] = flex.int(int(x[2]) + handle.starting_frame for x in xyzcal)

        # compute ZETA
        table.compute_zeta(self._experiment)
        # compute LP
        table.compute_corrections(self._experiments)

        table.set_flags(flex.bool(table.size(), True), table.flags.predicted)
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated)

        logger.info(f"Created table with {len(table)} reflections")

        # Output the table to pickle file
        if params.output.reflections is None:
            params.output.reflections = "integrate_hkl.refl"
        logger.info(f"Saving reflection table to {params.output.reflections}")
        table.as_file(params.output.reflections)
        logger.info(f"Saved reflection table to {params.output.reflections}")

    def derive_reindex_matrix(self, handle):
        """Derive a reindexing matrix to go from the orientation matrix used
        for XDS integration to the one used for DIALS integration."""
        dA = matrix.sqr(self._experiment.crystal.get_A())
        dbeam = matrix.col(self._experiment.beam.get_sample_to_source_direction())
        daxis = matrix.col(self._experiment.goniometer.get_rotation_axis())
        n = dbeam.cross(daxis)
        xbeam = matrix.col(handle.beam_vector).normalize()
        xaxis = matrix.col(handle.rotation_axis).normalize()

        # want to align XDS -s0 vector...
        R = align_reference_frame(-xbeam, dbeam, xaxis, n.cross(dbeam))
        xA = matrix.sqr(
            handle.unit_cell_a_axis + handle.unit_cell_b_axis + handle.unit_cell_c_axis
        ).inverse()
        xA = R * xA

        # assert that this should just be a simple integer rotation matrix
        # i.e. reassignment of a, b, c so...

        return matrix.sqr([int(round(e)) for e in (dA.inverse() * xA).elems])


class XDSFileImporter:
    """Import an experimentlist from xds.

    This will try to find the xds inp file - if that doesn't exist, it
    will try to find the best available file.
    """

    def __init__(self, xds_directory: Path):
        """Initialise with the options"""
        self.xds_directory = xds_directory

    def __call__(self, params, options):
        # Get the XDS.INP file
        xds_inp = self.xds_directory / "XDS.INP"
        if not xds_inp.exists():
            raise RuntimeError(f"Unable to find XDS.INP file in {self.xds_directory}.")

        if params.input.xds_file is None:
            xds_file = XDSFileImporter.find_best_xds_file(self.xds_directory)
        else:
            xds_file = params.input.xds_file

        # Check a file is given
        if xds_file is None:
            msg = "one of " + ", ".join(required_files_to_make_experiments)
            raise RuntimeError(f"No XDS file ({msg}) found in {self.xds_directory}")

        # Load the experiment list
        experiments = ExperimentListFactory.from_xds(xds_inp, xds_file)

        # set some dummy epochs
        nimg = (
            experiments[0].scan.get_image_range()[1]
            - experiments[0].scan.get_image_range()[0]
            + 1
        )
        epochs = flex.double(nimg, 0.0)
        for i in range(1, (nimg - 1)):
            epochs[i] = epochs[(i - 1)] + 1.0
        experiments[0].scan.set_epochs(epochs)
        experiments[0].scan.set_exposure_times(epochs)

        # Print some general info
        logger.info("-" * 80)
        logger.info(
            f"Read {len(experiments)} experiment{'s' if len(experiments) > 1 else ''} from {xds_file}"
        )

        # Attempt to create scan-varying crystal model if requested
        if params.read_varying_crystal:
            integrate_lp = os.path.join(self.args[0], "INTEGRATE.LP")
            if os.path.isfile(integrate_lp):
                self.extract_varying_crystal(integrate_lp, experiments)
            else:
                logger.info(
                    "No INTEGRATE.LP to extract varying crystal model. Skipping"
                )

        # Loop through the data blocks
        for i, exp in enumerate(experiments):
            # Print some experiment info
            logger.info("-" * 80)
            logger.info("Experiment %d" % i)
            logger.info(f"  format: {exp.imageset.get_format_class()}")
            logger.info(f"  type: {type(exp.imageset)}")
            logger.info(f"  num images: {len(exp.imageset)}")

            # Print some model info
            if options.verbose > 1:
                logger.info("")
                if exp.beam:
                    logger.info(exp.beam)
                else:
                    logger.info("no beam!")
                if exp.detector:
                    logger.info(exp.detector)
                else:
                    logger.info("no detector!")
                if exp.goniometer:
                    logger.info(exp.goniometer)
                else:
                    logger.info("no goniometer!")
                if exp.scan:
                    logger.info(exp.scan)
                else:
                    logger.info("no scan!")
                if exp.crystal:
                    logger.info(exp.crystal)
                else:
                    logger.info("no crystal!")

        # Write the experiment list
        logger.info("-" * 80)
        logger.info(f"Writing experiments to {params.output.xds_experiments}")
        experiments.as_file(params.output.xds_experiments)
        return experiments

    @staticmethod
    def find_best_xds_file(xds_dir):
        """Find the best available file."""

        # The possible files to check
        paths = [xds_dir / f for f in required_files_to_make_experiments]

        # Return the first path that exists
        for p in paths:
            if os.path.exists(p):
                return p

        # If no path exists, return None
        return None

    @staticmethod
    def extract_varying_crystal(integrate_lp, experiments):
        """Extract a varying crystal model from an INTEGRATE.LP file (static
        in blocks with step changes) and write it to the provided
        experiments
        """

        if len(experiments) > 1:
            logger.info(
                "Can only read a varying crystal model for a single "
                + "experiment. Skipping."
            )
            return
        experiment = experiments[0]

        # NOTE: we are not reading SEGMENT information here - so HKL
        # file needs also be without SEGMENT information, e.g. via
        #
        # egrep -v "^.SEGMENT=|^! .*=" INTEGRATE.HKL.orig | sed "s/,ISEG//g" | sed "s/^ \(.*\) [ ]*[0-9][0-9]*[ ]*$/ \1/g" > INTEGRATE.HKL

        # read required records from the file. Relies on them being in the
        # right order as we read through once
        xds_axis = None
        xds_beam = None
        blocks, a_axis, b_axis, c_axis = [], [], [], []
        with open(integrate_lp) as f:
            for record in f:
                if record.lstrip().startswith("ROTATION_AXIS="):
                    xds_axis = record.split("ROTATION_AXIS=")[1].split()
                    break
            for record in f:
                if record.lstrip().startswith("INCIDENT_BEAM_DIRECTION="):
                    xds_beam = record.split("INCIDENT_BEAM_DIRECTION=")[1].split()
                    break
            for record in f:
                if record.lstrip().startswith("PROCESSING OF IMAGES"):
                    blocks.append(record.split("PROCESSING OF IMAGES")[1])
                    continue
                if record.lstrip().startswith("COORDINATES OF UNIT CELL A-AXIS"):
                    a_axis.append(record.split("COORDINATES OF UNIT CELL A-AXIS")[1])
                    continue
                if record.lstrip().startswith("COORDINATES OF UNIT CELL B-AXIS"):
                    b_axis.append(record.split("COORDINATES OF UNIT CELL B-AXIS")[1])
                    continue
                if record.lstrip().startswith("COORDINATES OF UNIT CELL C-AXIS"):
                    c_axis.append(record.split("COORDINATES OF UNIT CELL C-AXIS")[1])
                    continue

        # sanity checks
        msg = "INTEGRATE.LP is not in the expected format"
        nblocks = len(blocks)
        try:
            assert len(a_axis) == len(b_axis) == len(c_axis) == nblocks
            assert (xds_axis, xds_beam).count(None) == 0
        except AssertionError:
            logger.info(msg)
            return

        # conversions to numeric
        try:
            blocks = [[int(b) for b in block.split("...")] for block in blocks]
            a_axis = [[float(a) for a in axis.split()] for axis in a_axis]
            b_axis = [[float(b) for b in axis.split()] for axis in b_axis]
            c_axis = [[float(c) for c in axis.split()] for axis in c_axis]
            xds_beam = [float(e) for e in xds_beam]
            xds_axis = [float(e) for e in xds_axis]
        except ValueError:
            logger.info(msg)
            return

        # coordinate frame conversions
        dbeam = matrix.col(experiment.beam.get_sample_to_source_direction())
        daxis = matrix.col(experiment.goniometer.get_rotation_axis())
        xbeam = matrix.col(xds_beam).normalize()
        xaxis = matrix.col(xds_axis).normalize()

        # want to align XDS -s0 vector...
        R = align_reference_frame(-xbeam, dbeam, xaxis, daxis)

        # Make a static crystal for each block
        crystals = []
        sg = experiment.crystal.get_space_group()
        for a, b, c in zip(a_axis, b_axis, c_axis):
            a = R * matrix.col(a)
            b = R * matrix.col(b)
            c = R * matrix.col(c)
            crystals.append(Crystal(a, b, c, space_group=sg))

        # construct a list of scan points
        A_list = []
        for block, crystal in zip(blocks, crystals):
            A = crystal.get_A()
            for im in range(block[0], block[1] + 1):
                A_list.append(A)
        # Need a final scan point at the end of the final image
        A_list.append(A)

        # set the scan-varying crystal
        experiment.crystal.set_A_at_scan_points(A_list)


phil_scope = parse("""
input {
    method = experiment reflections
        .type = choice
        .help = "Deprecated option - has no effect"

    xds_file = None
        .type = path
        .help = "Explicitly specify the file to use "
}

output {
    filename = None
        .type = str
        .help = "Deprecated option - has no effect"
    reflections = None
        .type = str
        .help = "The output filname of the reflections file (defaults to either integrate_hkl.refl or spot_xds.refl)"
    xds_experiments = "xds_models.expt"
        .type = str
        .help = "The output filename of the experiment list created from xds"
    log = dials.import_xds.log
        .type = path
}

remove_invalid = False
    .type = bool
    .help = "Remove non-index reflections (if miller indices are present)"

add_standard_columns = False
    .type = bool
    .help = "Add empty standard columns to the reflections. Note columns"
            "for centroid variances are set to contain 1s, not 0s"

read_varying_crystal = False
    .type = bool
    .help = "Attempt to create a scan-varying crystal model from"
            "INTEGRATE.LP, if present"
""")


@show_mail_handle_errors()
def run(args=None):
    usage = (
        "dials.import_xds /path/to/folder/containing/xds/inp/ (SPOT.XDS|INTEGRATE.HKL)"
    )

    parser = ArgumentParser(usage=usage, phil=phil_scope)
    # Parse the command line arguments
    params, options, unhandled = parser.parse_args(
        args, show_diff_phil=True, return_unhandled=True
    )

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.output.filename:
        logger.warning("""
The option output.filename= is deprecated and has no effect.
To set the output reflections filename, please use output.reflections=
To set the output experiment filename, please use output.xds_experiments=
""")

    if params.input.method:
        logger.warning("""
The option input.method= is deprecated and has no effect.
The program now attempts to output both experiment and reflection data
based on the input files available.
""")

    # Check number of arguments
    if len(unhandled) == 0:
        parser.print_help()
        exit(0)

    ## Assume that the unhandled arguments contain a directory and possibly a filepath or filename (e.g. "INTEGRATE.HKL", "SPOT.XDS")
    directories = [Path(u).resolve() for u in unhandled if Path(u).is_dir()]
    for d in list(directories):
        if d.name == "INTEGRATE.HKL" or d.name == "SPOT.XDS":
            directories.remove(d)

    if len(directories) == 0:
        # We might be in the implicit case where a path to the SPOT.XDS/INTEGRATE.HKL was provided (old usage)
        # if so, extract the directory
        directories = list({Path(u).resolve().parent for u in unhandled})

    if len(directories) != 1:
        msg = "\n ".join(str(d) for d in directories)
        logger.info(
            f"A single xds directory is required, found more than one in input arguments:\n {msg}"
        )
        sys.exit(0)

    xds_directory = Path(directories[0])

    # check whether we have an explicit choice of SPOT.XDS (to override default INTEGRATE.HKL)
    use_spot_xds: bool = False
    spot_xds: Optional[Path] = None
    integrate_hkl: Optional[Path] = None

    for arg in unhandled:
        if Path(arg).name == "SPOT.XDS":
            spot_xds = Path(arg)
        elif Path(arg).name == "INTEGRATE.HKL":
            integrate_hkl = Path(arg)

    if spot_xds and not integrate_hkl:
        if not spot_xds.exists():  # i.e. just the word on the command line
            if not (xds_directory / "SPOT.XDS").exists():
                logger.info(
                    f"Unable to find SPOT_XDS as specified or in {xds_directory}"
                )
                sys.exit(0)
            spot_xds = xds_directory / "SPOT.XDS"
        use_spot_xds = True

    if integrate_hkl:
        if not params.input.xds_file:  # use the specified integrate.hkl when creating the models, to ensure consistency.
            if integrate_hkl.is_file():
                params.input.xds_file = integrate_hkl
            elif (xds_directory / "INTEGRATE.HKL").is_file():
                params.input.xds_file = xds_directory / "INTEGRATE.HKL"
        if not integrate_hkl.exists():  # i.e. just the word on the command line
            if not (xds_directory / "INTEGRATE.HKL").exists():
                logger.info(
                    f"Unable to find INTEGRATE.HKL as specified or in {xds_directory}"
                )
                sys.exit(0)
            integrate_hkl = xds_directory / "INTEGRATE.HKL"

    # If neither explicitly specified on command line, default to using INTEGRATE.HKL, then SPOT.XDS if they exist
    if not spot_xds and not integrate_hkl:
        if (xds_directory / "INTEGRATE.HKL").exists():
            integrate_hkl = xds_directory / "INTEGRATE.HKL"
        elif (xds_directory / "SPOT.XDS").exists():
            spot_xds = xds_directory / "SPOT.XDS"
            use_spot_xds = True

    # First make the experiments
    importer = XDSFileImporter(xds_directory=xds_directory)
    try:
        expts = importer(params, options)
    except RuntimeError as e:
        logger.info(e)
        sys.exit(0)

    # If we specified SPOT.XDS (or only SPOT.XDS exists), use that
    if use_spot_xds:
        refl_importer = SpotXDSImporter(spot_xds)
    else:  # Else use the INTEGRATE.HKL as identified above.
        if not integrate_hkl:
            logger.info(
                "Unable to find SPOT.XDS or INTEGRATE.HKL in order to create a reflection table, finishing."
            )  # allow to just import expt file.
            sys.exit(0)
        refl_importer = IntegrateHKLImporter(integrate_hkl, expts)
    try:
        refl_importer(params, options)
    except RuntimeError as e:
        logger.info(e)
        sys.exit(0)


if __name__ == "__main__":
    run()
