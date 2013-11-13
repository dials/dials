#!/usr/bin/env python
#
# dials.algorithms.refinement.refiner2.py
#
#  Copyright (C) 2013 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Authors: James Parkhurst, David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Refiner is the refinement module public interface. RefinerFactory is what
should usually be used to construct a Refiner.

Here work on a new version of the interface"""

from __future__ import division
from dials.algorithms.refinement.refinement_helpers import print_model_geometry

class Bucket:
    pass

class RefinerFactory(object):
    """Factory class to create refiners"""

    @staticmethod
    def from_parameters_models_data(params,
                                    reflections,
                                    sweep=None,
                                    beam=None,
                                    goniometer=None,
                                    detector=None,
                                    scan=None,
                                    image_width_rad=None,
                                    sweep_range_rad=None,
                                    crystal=None,
                                    crystals=None,
                                    verbosity=0):
        """Given a set of parameters, experimental models and reflections,
        construct the refiner

        Params:
            params The input parameters
            verbosity The verbosity level

        Returns:
            The refiner instance

        NB it is possible to have image_width_rad=None, sweep_range_rad=None
        and scan=None and still do X, Y, Phi refinement, but in this case the
        PHIL parameters should not specify reflections.reflections_per_degree
        and must specify target.rmsd_cutoff=absolute.

        """

        # checks on the input
        if sweep:
            assert [beam, goniometer, detector, scan].count(None) == 4
            beam = sweep.get_beam()
            detector = sweep.get_detector()
            goniometer = sweep.get_goniometer()
            scan = sweep.get_scan()

        # if either of these provided, both must be
        if image_width_rad or sweep_range_rad:
            assert image_width_rad and sweep_range_rad

        # only one of image_width or scan should be provided (or neither)
        assert [image_width_rad, scan].count(None) >= 1
        if scan:
            image_width_rad = scan.get_oscillation(deg=False)[1]
            sweep_range_rad = scan.get_oscillation_range(deg=False)

        # do we have the essential models?
        assert [beam, detector].count(None) == 0
        assert [crystal, crystals].count(None) == 1

        # copy the models
        from dxtbx.model import Beam, Detector
        from dxtbx.array_family import flex
        import copy
        # use copy constructors
        beam = Beam(beam)
        detector = Detector(flex.panel([panel for panel in detector]))
        if crystal: crystals = [copy.deepcopy(crystal)]
        if crystals: crystals = copy.deepcopy(crystals)

        # copy the reflections
        reflections = reflections.deep_copy()

        # check that the beam vectors are stored: if not, compute them
        from scitbx import matrix
        for ref in reflections:
            if ref.beam_vector != (0.0, 0.0, 0.0):
                continue
            panel = detector[ref.panel_number]
            x, y = panel.millimeter_to_pixel(ref.image_coord_mm)
            ref.beam_vector = matrix.col(panel.get_pixel_lab_coord(
                (x, y))).normalize() / beam.get_wavelength()

        # create parameterisations
        parameterisations = \
                RefinerFactory.configure_parameterisation(
                    params, beam, detector, crystals, goniometer, scan)

        pred_param = parameterisations.pred_param

        if verbosity > 1:
            print ""
            print "Experimental Models"
            print "-------------------"
            print beam
            print detector
            if goniometer: print goniometer
            if scan: print scan
            for x in crystals: print x

        if verbosity > 1:
            print "Prediction equation parameterisation built\n"
            print "Parameter order : name mapping"
            for i, e in enumerate(pred_param.get_param_names()):
                print "Parameter %03d : " % i + e
            print

            print "Prior to refinement the experimental model is:"
            # FIXME. Won't work for a list of crystals. Do I really
            # need print_model_geometry? (probably not)
            print_model_geometry(beam, detector, crystals[0])
            print

        if verbosity > 1:
            print "Building reflection manager"
            print ("Input reflection list size = %d observations"
                   % len(reflections))

        # create reflection manager
        refman = RefinerFactory.configure_refman(params, reflections,
            beam, goniometer, sweep_range_rad, verbosity)

        if verbosity > 1:
            print ("Number of observations that pass inclusion criteria = %d"
                   % refman.get_total_size())
            print ("Working set size = %d observations"
                   % refman.get_sample_size())
            print "Reflection manager built\n"

        if verbosity > 1: print "Building target function"

        # create target function
        target = RefinerFactory.configure_target(params, crystals, beam,
                        goniometer, detector, image_width_rad, refman,
                        pred_param)

        if verbosity > 1: print "Target function built\n"

        if verbosity > 1: print "Building refinement engine"

        # create refinery
        refinery = RefinerFactory.configure_refinery(
                        params, target, pred_param, verbosity)

        if verbosity > 1: print "Refinement engine built\n"

        # build refiner interface and return
        return Refiner(beam, crystals, detector,
                       parameterisations, refman, target, refinery,
                       goniometer=goniometer,
                       scan=scan,
                       verbosity=verbosity)

    @staticmethod
    def configure_parameterisation(
            params, beam, detector, crystals, goniometer, scan):
        """Given a set of parameters, create a parameterisation from a set of
        experimental models.

        Params:
            params The input parameters

        Returns:
            The parameterisation instances in a bucket
        """

        # Shorten parameter paths
        beam_options = params.refinement.parameterisation.beam
        crystal_options = params.refinement.parameterisation.crystal
        detector_options = params.refinement.parameterisation.detector
        prediction_options = params.refinement.parameterisation.prediction

        # Shorten paths
        import dials.algorithms.refinement.parameterisation as par

        # Beam (accepts goniometer=None)
        beam_param = par.BeamParameterisationOrientation(beam, goniometer)
        if beam_options.fix:
            if beam_options.fix == "all":
                beam_param.set_fixed([True, True])
            elif beam_options.fix == "in_spindle_plane":
                beam_param.set_fixed([True, False])
            else: # can only get here if refinement.phil is broken
                raise RuntimeError("beam_options.fix value not recognised")

        # Crystal
        # FIXME:
        if len(crystals) > 1:
            raise RuntimeError("Multiple crystal parameterisation not"
                               "yet supported")
        crystal = crystals[0] # Reminder for FIXME
        if crystal_options.scan_varying:
            assert [goniometer, scan].count(None) == 0
            xl_ori_param = par.ScanVaryingCrystalOrientationParameterisation(
                                                crystal,
                                                scan.get_image_range(),
                                                crystal_options.num_intervals)
            xl_uc_param = par.ScanVaryingCrystalUnitCellParameterisation(
                                                crystal,
                                                scan.get_image_range(),
                                                crystal_options.num_intervals)
        else:
            xl_ori_param = par.CrystalOrientationParameterisation(crystal)
            xl_uc_param = par.CrystalUnitCellParameterisation(crystal)

        if crystal_options.fix:
            if crystal_options.fix == "all":
                xl_ori_param.set_fixed([True] * xl_ori_param.num_total())
                xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
            elif crystal_options.fix == "cell":
                xl_uc_param.set_fixed([True] * xl_uc_param.num_total())
            elif crystal_options.fix == "orientation":
                xl_ori_param.set_fixed([True] * xl_ori_param.num_total())
            else: # can only get here if refinement.phil is broken
                raise RuntimeError("crystal_options.fix value not recognised")

        # Detector
        if detector_options.panels == "automatic":
            if len(detector) > 1:
                det_param = par.DetectorParameterisationMultiPanel(detector, beam)
            else:
                det_param = par.DetectorParameterisationSinglePanel(detector)
        elif detector_options.panels == "single":
            det_param = DetectorParameterisationSinglePanel(detector)
        elif self._detector_par_options == "multiple":
            det_param = DetectorParameterisationMultiPanel(detector, beam)
        else: # can only get here if refinement.phil is broken
            raise RuntimeError("detector_options.panels value not recognised")

        if detector_options.fix:
            if detector_options.fix == "all":
                det_param.set_fixed([True] * det_param.num_total())
            elif detector_options.fix == "position":
                det_params = det_param.get_params(only_free = False)
                to_fix = [e.param_type.startswith('length') \
                          for e in det_params]
                det_param.set_fixed(to_fix)
            elif detector_options.fix == "orientation":
                det_params = det_param.get_params(only_free = False)
                to_fix = [e.param_type.startswith('angle') \
                          for e in det_params]
                det_param.set_fixed(to_fix)
            else: # can only get here if refinement.phil is broken
                raise RuntimeError("detector_options.fix value not recognised")

        # Prediction equation parameterisation
        crystal = crystals[0] # FIXME: multiple xls not yet supported
        if crystal_options.scan_varying:
            pred_param = par.VaryingCrystalPredictionParameterisation(
                detector, beam, crystal, goniometer,
                [det_param], [beam_param], [xl_ori_param], [xl_uc_param])
        elif goniometer is None:
            from dials.algorithms.refinement.single_shots.parameterisation import \
                XYPredictionParameterisation
            pred_param = XYPredictionParameterisation(
                detector, beam, crystal, goniometer,
                [det_param], [beam_param], [xl_ori_param], [xl_uc_param])
        else:
            pred_param = par.XYPhiPredictionParameterisation(
                detector, beam, crystal, goniometer,
                [det_param], [beam_param], [xl_ori_param], [xl_uc_param])

        # Parameter reporting
        param_reporter = par.ParameterReporter([det_param], [beam_param],
            [xl_ori_param], [xl_uc_param])

        parameterisations = Bucket()
        parameterisations.beam_param = beam_param
        parameterisations.xl_ori_param = xl_ori_param
        parameterisations.xl_uc_param = xl_uc_param
        parameterisations.det_param = det_param
        parameterisations.pred_param = pred_param
        parameterisations.param_reporter = param_reporter

        return parameterisations

    @staticmethod
    def configure_refinery(params, target, pred_param, verbosity):
        """Given a set of parameters, a target class, and a prediction
        parameterisation class, build a refinery

        Params:
            params The input parameters

        Returns:
            The refinery instance
        """

        # Shorten parameter path
        options = params.refinement.refinery

        import dials.algorithms.refinement.engine as engine

        if options.engine  == "SimpleLBFGS":
            from engine import SimpleLBFGS as refinery
        elif options.engine == "LBFGScurvs":
            from engine import LBFGScurvs as refinery
        elif options.engine == "GaussNewtonIterations":
            from engine import GaussNewtonIterations as refinery
        elif options.engine == "LevMarIterations":
            from engine import LevenbergMarquardtIterations as refinery
        else:
            raise RuntimeError("Refinement engine " + options.engine +
                               " not recognised")

        return refinery(target = target,
                        prediction_parameterisation = pred_param,
                        log = options.log,
                        verbosity = verbosity,
                        track_step = options.track_step,
                        track_gradient = options.track_gradient,
                        max_iterations = options.max_iterations)

    @staticmethod
    def configure_refman(params, reflections, beam, goniometer,
                         sweep_range_rad, verbosity):
        """Given a set of parameters and models, build a reflection manager

        Params:
            params The input parameters

        Returns:
            The reflection manager instance
        """

        # Shorten parameter path
        options = params.refinement.reflections
        options.random_seed
        if options.use_all_reflections:
            nref_per_degree = None
        else:
            nref_per_degree = options.reflections_per_degree

        # While a random subset of reflections is used, continue to
        # set random.seed to get consistent behaviour
        if options.random_seed:
            import random
            random.seed(options.random_seed)
            if verbosity > 1:
                print "Random seed set to %d\n" % options.random_seed

        if goniometer:
            import dials.algorithms.refinement.target as target
            from target import ReflectionManager as refman

        else:
            from dials.algorithms.refinement.single_shots.target import \
                ReflectionManagerXY as refman

        return refman(reflections=reflections,
                      beam=beam,
                      gonio=goniometer,
                      sweep_range_rad=sweep_range_rad,
                      nref_per_degree=nref_per_degree,
                      min_num_obs=options.minimum_number_of_reflections,
                      max_num_obs=options.maximum_number_of_reflections,
                      inclusion_cutoff=options.inclusion_cutoff,
                      verbosity=verbosity)

    @staticmethod
    def configure_target(params, crystals, beam, goniometer, detector,
        image_width_rad, refman, pred_param):
        """Given a set of parameters, configure a factory to build a
        target function

        Params:
            params The input parameters

        Returns:
            The target factory instance
        """

        # Shorten parameter path
        options = params.refinement.target

        if options.rmsd_cutoff == "fraction_of_bin_size":
            absolute_cutoffs = None
        elif options.rmsd_cutoff == "absolute":
            absolute_cutoffs = options.absolute_cutoffs
        else:
            raise RuntimeError("Target function rmsd_cutoff option" +
                options.rmsd_cutoff + " not recognised")

        # Determine whether the target is in X, Y, Phi space or just X, Y.
        crystal = crystals[0] # FIXME: multiple crystals not yet supported
        if goniometer:
            from dials.algorithms.refinement.prediction import ReflectionPredictor
            ref_predictor = ReflectionPredictor(crystal, beam, goniometer)

            import dials.algorithms.refinement.target as target
            from target import LeastSquaresPositionalResidualWithRmsdCutoff

            target = LeastSquaresPositionalResidualWithRmsdCutoff(
                            ref_predictor, detector, refman, pred_param,
                            image_width_rad, options.bin_size_fraction,
                            absolute_cutoffs)
        else:
            # FIXME: StillsReflectionPredictor doesn't actually
            # exist yet!
            from dials.algorithms.refinement.prediction import \
                StillsReflectionPredictor
            ref_predictor = ReflectionPredictor(crystal, beam)

            import dials.algorithms.refinement.single_shots.target as target
            from target import LeastSquaresXYResidualWithRmsdCutoff

            target = LeastSquaresXYResidualWithRmsdCutoff(
                            ref_predictor, detector, refman, pred_param,
                            options.bin_size_fraction,
                            absolute_cutoffs)

        return target

class Refiner(object):
    """The refiner class."""

    def __init__(self, beam, crystals, detector,
                 parameterisations, refman, target, refinery,
                 goniometer=None,
                 scan=None,
                 verbosity=0):
        """ Initialise the refiner class.

        Params:
            FIXME Documentation
        """

        # keep the models for access after refinement
        self.beam = beam
        self.crystals = crystals
        # only keep crystal if there is only one of them
        self.crystal = crystals[0] if len(crystals) == 1 else None
        self.goniometer = goniometer
        self.detector = detector

        # keep the parameterisations(FIXME should be private? - or
        # maybe don't even need these references to them?)
        self.beam_param = parameterisations.beam_param
        self.xl_ori_param = parameterisations.xl_ori_param
        self.xl_uc_param = parameterisations.xl_uc_param
        self.det_param = parameterisations.det_param

        # need this?
        self.pred_param = parameterisations.pred_param

        # I do want this though
        self.param_report = parameterisations.param_reporter

        # should be private?
        self.refman = refman
        self.target = target
        self.refinery = refinery

        self._verbosity = verbosity

        return

    def rmsds(self):
        """Return rmsds of the current model"""

        self.refinery.prepare_for_step()

        return self.target.rmsds()

    def __call__(self):
        """Run refinement"""

        ###################################
        # Do refinement and return models #
        ###################################

        self.refinery.run()

        if self._verbosity > 1:
            print
            print "Refinement has completed with the following geometry:"
            # FIXME Not useful for multiple crystals
            print_model_geometry(self.beam, self.detector, self.crystals[0])

            try: # can only do this if there is a scan
                if self.param_report.varying_params_vs_image_number(
                    self.scan.get_image_range()):
                    print "Writing scan-varying parameter table to file"
            except AttributeError:
                pass

            print "Reporting on the refined parameters:"
            print self.param_report

            print "Writing residuals to file"
            self.write_residuals_table()

        # Return the refinery, containing useful information such as the
        # refinement history. The refined models are set by side-effect
        return self.refinery

    def selection_used_for_refinement(self):
        """Return a selection as a flex.bool in terms of the input reflection
        data of those reflections that were used in the final step of
        refinement."""

        from scitbx.array_family import flex
        matches = self.refman.get_matches()
        selection = flex.bool(len(self.reflections))
        for m in matches:
            selection[m.iobs] = True

        return selection

    def write_residuals_table(self):

        matches = self.refman.get_matches()

        f = open("residuals.dat","w")
        header = ("H\tK\tL\tFrame_obs\tX_obs\tY_obs\tPhi_obs\tX_calc\t"
            "Y_calc\tPhi_calc\n")
        f.write(header)

        for m in matches:
            msg = ("%d\t%d\t%d\t%d\t%5.3f\t%5.3f\t%9.6f\t%5.3f\t%9.6f\t"
                  "%5.3f\n")
            msg = msg % (m.H[0], m.H[1], m.H[2], m.frame_o, m.Xo, m.Yo,
                         m.Phio, m.Xc, m.Yc, m.Phic)
            f.write(msg)
        f.close()

    def predict_reflections(self):
        """Predict all reflection positions after refinement"""

        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.spot_prediction import ray_intersection
        from dials.model.data import ReflectionList
        from math import pi
        from dials.algorithms.refinement.prediction.predictors import \
                ScanVaryingReflectionListGenerator

        s0 = self.beam.get_s0()
        dmin = self.detector.get_max_resolution(s0)
        sv_predictor = ScanVaryingReflectionListGenerator(self.pred_param,
                            self.beam, self.gonio, self.scan, dmin)

        # Duck typing to determine whether prediction is scan-varying or not
        try:
            refs = ReflectionList(sv_predictor())
            self._new_reflections = ray_intersection(self.detector, refs)

        except AttributeError: # prediction seems to be scan-static
            predict = ReflectionPredictor()
            self._new_reflections = predict(self.sweep, self.crystal)

        return self._new_reflections
