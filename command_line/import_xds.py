#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import os

from dials.array_family import flex


class SpotXDSImporter(object):
    """ Class to import a spot.xds file to a reflection table. """

    def __init__(self, spot_xds):
        self._spot_xds = spot_xds

    def __call__(self, params, options):
        """ Import the spot.xds file. """
        from iotbx.xds import spot_xds
        from dials.util.command_line import Command

        # Read the SPOT.XDS file
        Command.start("Reading SPOT.XDS")
        handle = spot_xds.reader()
        handle.read_file(self._spot_xds)
        centroid = handle.centroid
        intensity = handle.intensity
        try:
            miller_index = handle.miller_index
        except AttributeError:
            miller_index = None
        Command.end("Read {} spots from SPOT.XDS file.".format(len(centroid)))

        # Create the reflection list
        Command.start("Creating reflection list")
        table = flex.reflection_table()
        table["id"] = flex.int(len(centroid), 0)
        table["panel"] = flex.size_t(len(centroid), 0)
        if miller_index:
            table["miller_index"] = flex.miller_index(miller_index)
        table["xyzobs.px.value"] = flex.vec3_double(centroid)
        table["intensity.sum.value"] = flex.double(intensity)
        Command.end("Created reflection list")

        # Remove invalid reflections
        Command.start("Removing invalid reflections")
        if miller_index and params.remove_invalid:
            flags = flex.bool([h != (0, 0, 0) for h in table["miller_index"]])
            table = table.select(flags)
        Command.end("Removed invalid reflections, %d remaining" % len(table))

        # Fill empty standard columns
        if params.add_standard_columns:
            Command.start("Adding standard columns")
            rt = flex.reflection_table.empty_standard(len(table))
            rt.update(table)
            table = rt
            # set variances to unity
            table["xyzobs.mm.variance"] = flex.vec3_double(len(table), (1, 1, 1))
            table["xyzobs.px.variance"] = flex.vec3_double(len(table), (1, 1, 1))
            Command.end("Standard columns added")

        # Output the table to pickle file
        if params.output.filename is None:
            params.output.filename = "spot_xds.refl"
        Command.start("Saving reflection table to %s" % params.output.filename)
        table.as_pickle(params.output.filename)
        Command.end("Saved reflection table to %s" % params.output.filename)


class IntegrateHKLImporter(object):
    """ Class to import an integrate.hkl file to a reflection table. """

    def __init__(self, integrate_hkl, experiment):
        self._integrate_hkl = integrate_hkl
        self._experiment = experiment

    def __call__(self, params, options):
        """ Import the integrate.hkl file. """

        from iotbx.xds import integrate_hkl
        from dials.util.command_line import Command
        from cctbx import sgtbx

        # Get the unit cell to calculate the resolution
        uc = self._experiment.crystal.get_unit_cell()

        # Read the INTEGRATE.HKL file
        Command.start("Reading INTEGRATE.HKL")
        handle = integrate_hkl.reader()
        handle.read_file(self._integrate_hkl)
        hkl = flex.miller_index(handle.hkl)
        xyzcal = flex.vec3_double(handle.xyzcal)
        xyzobs = flex.vec3_double(handle.xyzobs)
        iobs = flex.double(handle.iobs)
        sigma = flex.double(handle.sigma)
        rlp = flex.double(handle.rlp)
        peak = flex.double(handle.peak) * 0.01
        if len(handle.iseg):
            panel = flex.size_t(handle.iseg) - 1
        else:
            panel = flex.size_t(len(hkl), 0)
        Command.end("Read %d reflections from INTEGRATE.HKL file." % len(hkl))

        if len(self._experiment.detector) > 1:
            for p_id, p in enumerate(self._experiment.detector):
                sel = panel == p_id
                offset = p.get_raw_image_offset()
                xyzcal.set_selected(sel, xyzcal.select(sel) - (offset[0], offset[1], 0))
                xyzobs.set_selected(sel, xyzobs.select(sel) - (offset[0], offset[1], 0))

        # Derive the reindex matrix
        rdx = self.derive_reindex_matrix(handle)
        print("Reindex matrix:\n%d %d %d\n%d %d %d\n%d %d %d" % (rdx.elems))

        # Reindex the reflections
        Command.start("Reindexing reflections")
        cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(sgtbx.rot_mx(rdx.elems)))
        hkl = cb_op.apply(hkl)
        Command.end("Reindexed %d reflections" % len(hkl))

        # Create the reflection list
        Command.start("Creating reflection table")
        table = flex.reflection_table()
        table["id"] = flex.int(len(hkl), 0)
        table["panel"] = panel
        table["miller_index"] = hkl
        table["xyzcal.px"] = xyzcal
        table["xyzobs.px.value"] = xyzobs
        table["intensity.cor.value"] = iobs
        table["intensity.cor.variance"] = sigma ** 2
        table["intensity.prf.value"] = iobs * peak / rlp
        table["intensity.prf.variance"] = (sigma * peak / rlp) ** 2
        table["lp"] = 1.0 / rlp
        table["d"] = flex.double(uc.d(h) for h in hkl)
        Command.end("Created table with {} reflections".format(len(table)))

        # Output the table to pickle file
        if params.output.filename is None:
            params.output.filename = "integrate_hkl.refl"
        Command.start("Saving reflection table to %s" % params.output.filename)
        table.as_pickle(params.output.filename)
        Command.end("Saved reflection table to %s" % params.output.filename)

    def derive_reindex_matrix(self, handle):
        """Derive a reindexing matrix to go from the orientation matrix used
        for XDS integration to the one used for DIALS integration."""
        from scitbx import matrix

        dA = matrix.sqr(self._experiment.crystal.get_A())
        dbeam = matrix.col(self._experiment.beam.get_direction())
        daxis = matrix.col(self._experiment.goniometer.get_rotation_axis())
        n = dbeam.cross(daxis)
        xbeam = matrix.col(handle.beam_vector).normalize()
        xaxis = matrix.col(handle.rotation_axis).normalize()

        # want to align XDS -s0 vector...
        from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

        R = align_reference_frame(-xbeam, dbeam, xaxis, n.cross(dbeam))
        xA = matrix.sqr(
            handle.unit_cell_a_axis + handle.unit_cell_b_axis + handle.unit_cell_c_axis
        ).inverse()
        xA = R * xA

        # assert that this should just be a simple integer rotation matrix
        # i.e. reassignment of a, b, c so...

        return matrix.sqr([int(round(e)) for e in (dA.inverse() * xA).elems])


class XDSFileImporter(object):
    """ Import a data block from xds. """

    def __init__(self, args):
        """ Initialise with the options"""
        self.args = args

    def __call__(self, params, options):
        from dxtbx.model.experiment_list import ExperimentListFactory
        from dxtbx.model.experiment_list import ExperimentListDumper

        # Get the XDS.INP file
        xds_inp = os.path.join(self.args[0], "XDS.INP")
        if params.input.xds_file is None:
            xds_file = XDSFileImporter.find_best_xds_file(self.args[0])
        else:
            xds_file = os.path.join(self.args[0], params.input.xds_file)

        # Check a file is given
        if xds_file is None:
            raise RuntimeError("No XDS file found")

        # Load the experiment list
        unhandled = []
        experiments = ExperimentListFactory.from_xds(xds_inp, xds_file)

        # Print out any unhandled files
        if len(unhandled) > 0:
            print("-" * 80)
            print("The following command line arguments were not handled:")
            for filename in unhandled:
                print("  %s" % filename)

        # Print some general info
        print("-" * 80)
        print("Read %d experiments from %s" % (len(experiments), xds_file))

        # Attempt to create scan-varying crystal model if requested
        if params.read_varying_crystal:
            integrate_lp = os.path.join(self.args[0], "INTEGRATE.LP")
            if os.path.isfile(integrate_lp):
                self.extract_varying_crystal(integrate_lp, experiments)
            else:
                print("No INTEGRATE.LP to extract varying crystal model. Skipping")

        # Loop through the data blocks
        for i, exp in enumerate(experiments):

            # Print some experiment info
            print("-" * 80)
            print("Experiment %d" % i)
            print("  format: %s" % str(exp.imageset.get_format_class()))
            print("  type: %s" % type(exp.imageset))
            print("  num images: %d" % len(exp.imageset))

            # Print some model info
            if options.verbose > 1:
                print("")
                if exp.beam:
                    print(exp.beam)
                else:
                    print("no beam!")
                if exp.detector:
                    print(exp.detector)
                else:
                    print("no detector!")
                if exp.goniometer:
                    print(exp.goniometer)
                else:
                    print("no goniometer!")
                if exp.scan:
                    print(exp.scan)
                else:
                    print("no scan!")
                if exp.crystal:
                    print(exp.crystal)
                else:
                    print("no crystal!")

        # Write the experiment list to a JSON or pickle file
        if params.output.filename is None:
            params.output.filename = "xds_models.expt"
        print("-" * 80)
        print("Writing experiments to %s" % params.output.filename)
        dump = ExperimentListDumper(experiments)
        dump.as_file(params.output.filename)

        # Optionally save as a data block
        if params.output.xds_experiments:
            print("-" * 80)
            print("Writing data block to %s" % params.output.xds_experiments)
            dump = ExperimentListDumper(experiments)
            dump.as_file(params.output.xds_experiments)

    @staticmethod
    def find_best_xds_file(xds_dir):
        """ Find the best available file."""

        # The possible files to check
        paths = [
            os.path.join(xds_dir, "XDS_ASCII.HKL"),
            os.path.join(xds_dir, "INTEGRATE.HKL"),
            os.path.join(xds_dir, "GXPARM.XDS"),
            os.path.join(xds_dir, "XPARM.XDS"),
        ]

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
            print(
                "Can only read a varying crystal model for a single "
                + "experiment. Skipping."
            )
            return
        experiment = experiments[0]

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
            print(msg)
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
            print(msg)
            return

        # coordinate frame conversions
        from scitbx import matrix

        dbeam = matrix.col(experiment.beam.get_direction())
        daxis = matrix.col(experiment.goniometer.get_rotation_axis())
        xbeam = matrix.col(xds_beam).normalize()
        xaxis = matrix.col(xds_axis).normalize()

        # want to align XDS -s0 vector...
        from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

        R = align_reference_frame(-xbeam, dbeam, xaxis, daxis)

        # Make a static crystal for each block
        from dxtbx.model import Crystal

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


class Script(object):
    """ A class to encapsulate the script. """

    def __init__(self):
        """ Initialise the script. """
        from dials.util.options import OptionParser
        from libtbx.phil import parse
        import libtbx.load_env

        # Create the phil parameters
        phil_scope = parse(
            """
      input {
        method = *experiment reflections
          .type = choice
          .help = "The input method"

        xds_file = None
          .type = str
          .help = "Explicitly specify the file to use"
      }

      output {
        filename = None
          .type = str
          .help = "The output file"

        xds_experiments = None
          .type = str
          .help = "Output filename of data block with xds"
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
    """
        )

        # The option parser
        usage = (
            "usage: %s [options] (SPOT.XDS|INTEGRATE.HKL)" % libtbx.env.dispatcher_name
        )
        self.parser = OptionParser(usage=usage, phil=phil_scope)

    def run(self):
        """ Run the script. """

        # Parse the command line arguments
        params, options, args = self.parser.parse_args(
            show_diff_phil=True, return_unhandled=True
        )

        # Check number of arguments
        if len(args) == 0:
            self.parser.print_help()
            exit(0)

        # Select the importer class
        if params.input.method == "experiment":
            importer = XDSFileImporter(args)
        else:
            importer = self.select_importer(args)

        # Import the XDS data
        importer(params, options)

    def select_importer(self, args):
        from dxtbx.model.experiment_list import ExperimentListFactory

        path, filename = os.path.split(args[0])
        if filename == "SPOT.XDS":
            return SpotXDSImporter(args[0])
        elif filename == "INTEGRATE.HKL":
            assert len(args) == 2
            experiments = ExperimentListFactory.from_json_file(args[1])
            assert len(experiments) == 1
            return IntegrateHKLImporter(args[0], experiments[0])
        else:
            raise RuntimeError("expected (SPOT.XDS|INTEGRATE.HKL), got %s" % filename)


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
