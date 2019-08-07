#!/usr/bin/env python
#
#  two_theta_refine.py
#
#  Copyright (C) 2016 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import datetime
import logging
import math
from time import time

from dials.array_family import flex
from dials.util import log
from dials.util.version import dials_version
from libtbx.utils import format_float_with_standard_uncertainty
from dials.util import Sorry
from dials.algorithms.refinement.corrgram import create_correlation_plots

logger = logging.getLogger("dials.command_line.two_theta_refine")

help_message = """

Refine the unit cell(s) of input experiments against the input indexed
reflections using a 2theta angle target. Report the refined cell and its
estimated standard deviation.

Examples::

  dials.two_theta_refine integrated.expt integrated.refl

  dials.two_theta_refine integrated.expt integrated.refl \
    correlation_plot.filename=corrplot.png cif=refined_cell.cif
"""

# The phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {
    experiments = refined_cell.expt
      .type = str
      .help = "The filename for experimental models including refined cells"

    log = dials.two_theta_refine.log
      .type = str

    debug_log = dials.two_theta_refine.debug.log
      .type = str

    cif = None
      .type = str
      .help = "Write unit cell error information to a Crystallographic"
              "Information File (CIF)"

    mmcif = None
      .type = str
      .help = "Write unit cell error information to a macromolecular"
              "Crystallographic Information File (mmCIF)"

    p4p = None
      .type = str
      .help = "Write output to SHELX / XPREP .p4p file"

    include scope dials.algorithms.refinement.corrgram.phil_scope
  }

  #FIXME expose _some_ of the Refiner options?
  #include scope dials.algorithms.refinement.refiner.phil_scope

  refinement
    .help = "Parameters to configure the refinement"
  {
    filter_integrated_centroids = True
      .type = bool
      .help = "If integrated centroids are provided, filter these so that only"
              "those with both the 'integrated' and 'strong' flags are used"

    combine_crystal_models = True
      .type = bool
      .help = "When multiple experiments are provided as input, combine these to"
              "fit the best single crystal model for all the data, or keep these"
              "models separate."

    triclinic = False
      .type = bool
      .help = "If true remove symmetry constraints and refine a triclinic cell"
              "by converting to P 1"
  }
""",
    process_includes=True,
)

working_phil = phil_scope.fetch()


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "models.expt observations.refl" % libtbx.env.dispatcher_name
        )

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=working_phil,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )

    @staticmethod
    def check_input(reflections):
        """Check the input is suitable for refinement. So far just check keys in
        the reflection table. Maybe later check experiments have overlapping models
        etc."""

        msg = (
            "The supplied reflection table does not have the required data "
            + "column: {0}"
        )
        for key in ["xyzobs.mm.value", "xyzobs.mm.variance"]:
            if key not in reflections:
                msg = msg.format(key)
                raise Sorry(msg)

        # FIXME add other things to be checked here
        return

    @staticmethod
    def combine_crystals(experiments):
        """Replace all crystals in the experiments list with the first crystal"""

        from dxtbx.model.experiment_list import Experiment, ExperimentList

        new_experiments = ExperimentList()
        ref_crystal = experiments[0].crystal
        for exp in experiments:
            new_experiments.append(
                Experiment(
                    beam=exp.beam,
                    detector=exp.detector,
                    scan=exp.scan,
                    goniometer=exp.goniometer,
                    crystal=ref_crystal,
                    imageset=exp.imageset,
                    identifier=exp.identifier,
                )
            )
        return new_experiments

    @staticmethod
    def filter_integrated_centroids(reflections):
        """Filter reflections to include only those with the integrated and the
        strong flag set, but only if there are apparently some integrated
        reflections"""

        orig_len = len(reflections)
        mask = reflections.get_flags(reflections.flags.integrated)
        if mask.count(True) == 0:
            return reflections
        reflections = reflections.select(mask)
        mask = reflections.get_flags(reflections.flags.strong)
        reflections = reflections.select(mask)

        logger.info(
            "{0} out of {1} reflections remain after filtering to keep only strong"
            " and integrated centroids".format(len(reflections), orig_len)
        )
        return reflections

    @staticmethod
    def convert_to_P1(reflections, experiments):
        """Convert the input crystals to P 1 and reindex the reflections"""
        from cctbx.sgtbx import space_group

        for iexp, exp in enumerate(experiments):
            sel = reflections["id"] == iexp
            xl = exp.crystal
            sg = xl.get_space_group()
            op = sg.info().change_of_basis_op_to_primitive_setting()
            exp.crystal = xl.change_basis(op)
            exp.crystal.set_space_group(space_group("P 1"))
            hkl_reindexed = op.apply(reflections["miller_index"].select(sel))
            reflections["miller_index"].set_selected(sel, hkl_reindexed)
        return reflections, experiments

    @staticmethod
    def create_refiner(params, reflections, experiments):

        from dials.algorithms.refinement.parameterisation.crystal_parameters import (
            CrystalUnitCellParameterisation,
        )
        from dials.algorithms.refinement.parameterisation.parameter_report import (
            ParameterReporter,
        )
        from dials.algorithms.refinement.two_theta_refiner import (
            TwoThetaReflectionManager,
            TwoThetaTarget,
            TwoThetaExperimentsPredictor,
            TwoThetaPredictionParameterisation,
        )

        # Only parameterise the crystal unit cell
        det_params = None
        beam_params = None
        xlo_params = None
        xluc_params = []
        for icrystal, crystal in enumerate(experiments.crystals()):
            exp_ids = experiments.indices(crystal)
            xluc_params.append(
                CrystalUnitCellParameterisation(crystal, experiment_ids=exp_ids)
            )

        # Two theta prediction equation parameterisation
        pred_param = TwoThetaPredictionParameterisation(
            experiments, det_params, beam_params, xlo_params, xluc_params
        )
        param_reporter = ParameterReporter(
            det_params, beam_params, xlo_params, xluc_params
        )

        # ReflectionManager, currently without outlier rejection
        # Note: If not all reflections are used, then the filtering must be
        # communicated to generate_cif/mmcif() to be included in the CIF file!
        refman = TwoThetaReflectionManager(
            reflections, experiments, outlier_detector=None
        )

        # Reflection predictor
        ref_predictor = TwoThetaExperimentsPredictor(experiments)

        # Two theta target
        target = TwoThetaTarget(experiments, ref_predictor, refman, pred_param)

        # Switch on correlation matrix tracking if a correlation plot is requested
        journal = None
        if params.output.correlation_plot.filename is not None:
            from dials.algorithms.refinement.engine import refinery_phil_scope

            journal = refinery_phil_scope.extract().refinery.journal
            journal.track_parameter_correlation = True

        # Minimisation engine - hardcoded to LevMar for now.
        from dials.algorithms.refinement.engine import (
            LevenbergMarquardtIterations as Refinery,
        )

        refinery = Refinery(
            target=target,
            prediction_parameterisation=pred_param,
            log=None,
            tracking=journal,
            max_iterations=20,
        )

        # Refiner
        from dials.algorithms.refinement.refiner import Refiner

        refiner = Refiner(
            experiments=experiments,
            pred_param=pred_param,
            param_reporter=param_reporter,
            refman=refman,
            target=target,
            refinery=refinery,
        )

        return refiner

    @staticmethod
    def cell_param_table(crystal):
        """Construct a table of cell parameters and their ESDs"""

        from libtbx.table_utils import simple_table

        cell = crystal.get_unit_cell().parameters()
        esd = crystal.get_cell_parameter_sd()
        vol = crystal.get_unit_cell().volume()
        vol_esd = crystal.get_cell_volume_sd()
        header = ["Parameter", "Value", "Estimated sd"]
        rows = []
        names = ["a", "b", "c", "alpha", "beta", "gamma"]
        for n, p, e in zip(names, cell, esd):
            rows.append([n, "%9.5f" % p, "%9.5f" % e])
        rows.append(["\nvolume", "\n%9.5f" % vol, "\n%9.5f" % vol_esd])
        st = simple_table(rows, header)
        return st.format()

    @staticmethod
    def generate_p4p(crystal, beam, file):
        logger.info("Saving P4P info to %s" % file)
        cell = crystal.get_unit_cell().parameters()
        esd = crystal.get_cell_parameter_sd()
        vol = crystal.get_unit_cell().volume()
        vol_esd = crystal.get_cell_volume_sd()

        open(file, "w").write(
            "\n".join(
                [
                    "TITLE    Auto-generated .p4p file from dials.two_theta_refine",
                    "CELL     %.4f %.4f %.4f %.4f %.4f %.4f %.4f"
                    % tuple(cell + (vol,)),
                    "CELLSD   %.4f %.4f %.4f %.4f %.4f %.4f %.4f"
                    % tuple(esd + (vol_esd,)),
                    "SOURCE   SYNCH   %.6f" % beam.get_wavelength(),
                    "",
                ]
            )
        )
        return

    @staticmethod
    def generate_cif(crystal, refiner, file):
        logger.info("Saving CIF information to %s" % file)
        from cctbx import miller
        import iotbx.cif.model

        block = iotbx.cif.model.block()
        block["_audit_creation_method"] = dials_version()
        block["_audit_creation_date"] = datetime.date.today().isoformat()
        #   block["_publ_section_references"] = '' # once there is a reference...

        for cell, esd, cifname in zip(
            crystal.get_unit_cell().parameters(),
            crystal.get_cell_parameter_sd(),
            [
                "length_a",
                "length_b",
                "length_c",
                "angle_alpha",
                "angle_beta",
                "angle_gamma",
            ],
        ):
            block["_cell_%s" % cifname] = format_float_with_standard_uncertainty(
                cell, esd
            )
        block["_cell_volume"] = format_float_with_standard_uncertainty(
            crystal.get_unit_cell().volume(), crystal.get_cell_volume_sd()
        )

        used_reflections = refiner.get_matches()
        block["_cell_measurement_reflns_used"] = len(used_reflections)
        block["_cell_measurement_theta_min"] = (
            flex.min(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_cell_measurement_theta_max"] = (
            flex.max(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_diffrn_reflns_number"] = len(used_reflections)
        miller_span = miller.index_span(used_reflections["miller_index"])
        min_h, min_k, min_l = miller_span.min()
        max_h, max_k, max_l = miller_span.max()
        block["_diffrn_reflns_limit_h_min"] = min_h
        block["_diffrn_reflns_limit_h_max"] = max_h
        block["_diffrn_reflns_limit_k_min"] = min_k
        block["_diffrn_reflns_limit_k_max"] = max_k
        block["_diffrn_reflns_limit_l_min"] = min_l
        block["_diffrn_reflns_limit_l_max"] = max_l
        block["_diffrn_reflns_theta_min"] = (
            flex.min(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_diffrn_reflns_theta_max"] = (
            flex.max(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )

        cif = iotbx.cif.model.cif()
        cif["two_theta_refine"] = block
        with open(file, "w") as fh:
            cif.show(out=fh)

    @staticmethod
    def generate_mmcif(crystal, refiner, file):
        logger.info("Saving mmCIF information to %s" % file)
        from cctbx import miller
        import iotbx.cif.model

        block = iotbx.cif.model.block()
        block["_audit.creation_method"] = dials_version()
        block["_audit.creation_date"] = datetime.date.today().isoformat()
        #   block["_publ.section_references"] = '' # once there is a reference...

        for cell, esd, cifname in zip(
            crystal.get_unit_cell().parameters(),
            crystal.get_cell_parameter_sd(),
            [
                "length_a",
                "length_b",
                "length_c",
                "angle_alpha",
                "angle_beta",
                "angle_gamma",
            ],
        ):
            block["_cell.%s" % cifname] = "%.8f" % cell
            block["_cell.%s_esd" % cifname] = "%.8f" % esd
        block["_cell.volume"] = "%f" % crystal.get_unit_cell().volume()
        block["_cell.volume_esd"] = "%f" % crystal.get_cell_volume_sd()

        used_reflections = refiner.get_matches()
        block["_cell_measurement.reflns_used"] = len(used_reflections)
        block["_cell_measurement.theta_min"] = (
            flex.min(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_cell_measurement.theta_max"] = (
            flex.max(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_diffrn_reflns.number"] = len(used_reflections)
        miller_span = miller.index_span(used_reflections["miller_index"])
        min_h, min_k, min_l = miller_span.min()
        max_h, max_k, max_l = miller_span.max()
        block["_diffrn_reflns.limit_h_min"] = min_h
        block["_diffrn_reflns.limit_h_max"] = max_h
        block["_diffrn_reflns.limit_k_min"] = min_k
        block["_diffrn_reflns.limit_k_max"] = max_k
        block["_diffrn_reflns.limit_l_min"] = min_l
        block["_diffrn_reflns.limit_l_max"] = max_l
        block["_diffrn_reflns.theta_min"] = (
            flex.min(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )
        block["_diffrn_reflns.theta_max"] = (
            flex.max(used_reflections["2theta_obs.rad"]) * 180 / math.pi / 2
        )

        cif = iotbx.cif.model.cif()
        cif["two_theta_refine"] = block
        with open(file, "w") as fh:
            cif.show(out=fh)

    def run(self):
        """Execute the script."""
        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)

        # set up global experiments and reflections lists
        from dials.array_family import flex

        reflections = flex.reflection_table()
        global_id = 0
        from dxtbx.model.experiment_list import ExperimentList

        experiments = ExperimentList()

        # loop through the input, building up the global lists
        nrefs_per_exp = []
        for ref_wrapper, exp_wrapper in zip(
            params.input.reflections, params.input.experiments
        ):
            refs = ref_wrapper.data
            exps = exp_wrapper.data
            for i, exp in enumerate(exps):
                sel = refs["id"] == i
                sub_ref = refs.select(sel)
                nrefs_per_exp.append(len(sub_ref))
                sub_ref["id"] = flex.int(len(sub_ref), global_id)
                reflections.extend(sub_ref)
                experiments.append(exp)
                global_id += 1

        # Try to load the models and data
        nexp = len(experiments)
        if nexp == 0:
            print("No Experiments found in the input")
            self.parser.print_help()
            return
        if len(reflections) == 0:
            print("No reflection data found in the input")
            self.parser.print_help()
            return

        self.check_input(reflections)

        # Configure the logging
        log.config(info=params.output.log, debug=params.output.debug_log)
        logger.info(dials_version())

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Convert to P 1?
        if params.refinement.triclinic:
            reflections, experiments = self.convert_to_P1(reflections, experiments)

        # Combine crystals?
        if params.refinement.combine_crystal_models and len(experiments) > 1:
            logger.info("Combining {0} crystal models".format(len(experiments)))
            experiments = self.combine_crystals(experiments)

        # Filter integrated centroids?
        if params.refinement.filter_integrated_centroids:
            reflections = self.filter_integrated_centroids(reflections)

        # Get the refiner
        logger.info("Configuring refiner")
        refiner = self.create_refiner(params, reflections, experiments)

        # Refine the geometry
        if nexp == 1:
            logger.info("Performing refinement of a single Experiment...")
        else:
            logger.info("Performing refinement of {} Experiments...".format(nexp))
        refiner.run()

        # get the refined experiments
        experiments = refiner.get_experiments()
        crystals = experiments.crystals()

        if len(crystals) == 1:
            # output the refined model for information
            logger.info("")
            logger.info("Final refined crystal model:")
            logger.info(crystals[0])
            logger.info(self.cell_param_table(crystals[0]))

        # Save the refined experiments to file
        output_experiments_filename = params.output.experiments
        logger.info(
            "Saving refined experiments to {}".format(output_experiments_filename)
        )
        from dxtbx.model.experiment_list import ExperimentListDumper

        dump = ExperimentListDumper(experiments)
        dump.as_json(output_experiments_filename)

        # Create correlation plots
        if params.output.correlation_plot.filename is not None:
            create_correlation_plots(refiner, params.output)

        if params.output.cif is not None:
            self.generate_cif(crystals[0], refiner, file=params.output.cif)

        if params.output.p4p is not None:
            self.generate_p4p(crystals[0], experiments[0].beam, file=params.output.p4p)

        if params.output.mmcif is not None:
            self.generate_mmcif(crystals[0], refiner, file=params.output.mmcif)

        # Log the total time taken
        logger.info("\nTotal time taken: {:.2f}s".format(time() - start_time))


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
