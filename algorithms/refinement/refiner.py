#!/usr/bin/env python
#
# dials.algorithms.refinement.refiner.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Authors: James Parkhurst, David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.algorithms.refinement import print_model_geometry


class Refiner(object):
    ''' The refiner class. '''

    def __init__(self, parameterisation_strategy, refinery_strategy,
                 reflections_strategy, target_strategy, verbosity):
        ''' Initialise the refiner class.

        Params:
            parameterisation_strategy The parameterisation strategy
            refinery_strategy The engine strategy
            reflections_strategy The reflection manager strategy
            target_strategy The target function strategy
            verbosity The verbosity level
        '''

        self.create_param = parameterisation_strategy
        self.create_refinery = refinery_strategy
        self.create_target = target_strategy
        self.create_refman = reflections_strategy
        self._verbosity = verbosity

    def prepare(self, sweep, crystal, reflections):
        ''' Prepare refiner with experimental models and data.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflection list

        Returns:
            The refined (sweep, crystal)

        '''
        from scitbx import matrix
        from math import sqrt
        import random

        # While a random subset of reflections is used, continue to
        # set random.seed to get consistent behaviour
        random.seed(42)
        if self._verbosity > 1: print "Random seed set to 42\n"

        # Get the models from the sweep
        self.beam = sweep.get_beam()
        self.detector = sweep.get_detector()
        self.gonio = sweep.get_goniometer()
        self.scan = sweep.get_scan()
        self.crystal = crystal

        if self._verbosity > 1:
            print ""
            print "Experimental Models"
            print "-------------------"
            print self.beam
            print self.detector
            print self.gonio
            print self.scan
            print self.crystal


        # Copy the reflections
        self.reflections = reflections
        self._saved_reflections = self.reflections.deep_copy()

        # check that the beam vectors are stored: if not, compute them
        # perhaps we should just be computing them here...
        for ref in self.reflections:
            if ref.beam_vector != (0.0, 0.0, 0.0):
                continue
            x, y = self.detector.millimeter_to_pixel(ref.image_coord_mm)
            ref.beam_vector = matrix.col(self.detector.get_pixel_lab_coord(
                (x, y))).normalize() / self.beam.get_wavelength()

        ###########################
        # Parameterise the models #
        ###########################

        (self.beam_param, self.xl_ori_param, self.xl_uc_param, self.det_param,
         self.pred_param) = \
            self.create_param(self.beam, self.crystal, self.gonio,
                              self.detector, self.scan)

        if self._verbosity > 1:
            print "Prediction equation parameterisation built\n"
            print "Parameter order : name mapping"
            for i, e in enumerate(self.pred_param.get_param_names()):
                print "Parameter %03d : " % i + e
            print

            print "Prior to refinement the experimental model is:"
            print_model_geometry(self.beam, self.detector, self.crystal)
            print

        #####################################
        # Select reflections for refinement #
        #####################################

        if self._verbosity > 1:
            print "Building reflection manager"
            print "Input reflection list size = %d observations" % len(self.reflections)

        self.refman = self.create_refman(self.reflections, self.beam,
                                self.gonio, self.scan, self._verbosity)

        if self._verbosity > 1:
            print "Number of observations that pass inclusion criteria = %d" % self.refman.get_total_size()
            print "Working set size = %d observations" % self.refman.get_sample_size()
            print "Reflection manager built\n"

        ##############################
        # Set up the target function #
        ##############################

        if self._verbosity > 1: print "Building target function"

        self.target = self.create_target(self.crystal, self.beam,
            self.gonio, self.detector, self.scan, self.refman,
            self.pred_param)

        if self._verbosity > 1: print "Target function built\n"

        ################################
        # Set up the refinement engine #
        ################################

        if self._verbosity > 1: print "Building refinement engine"

        self.refinery = self.create_refinery(self.target, self.pred_param,
                                             self._verbosity)

        if self._verbosity > 1: print "Refinement engine built\n"

        return

    def rmsds(self):
        '''Return rmsds of the current model'''

        self.refinery.prepare_for_step()

        return self.target.rmsds()

    def __call__(self, sweep=None, crystal=None, reflections=None):
        '''Run refinement'''

        if sweep and crystal and reflections:
            self.prepare(sweep, crystal, reflections)
        else: assert [sweep, crystal, reflections].count(None) == 3

        ###################################
        # Do refinement and return models #
        ###################################

        self.refinery.run()

        if self._verbosity > 1:
            print
            print "Refinement has completed with the following geometry:"
            print_model_geometry(self.beam, self.detector, self.crystal)

        # Do a test of new reflection pos
        #self._update_reflections_test()

        # Return the refinery, containing useful information such as the
        # refinement history. The refined models are set by side-effect
        return self.refinery

    def _update_reflections_test(self):
        from cctbx.array_family import flex
        from collections import defaultdict

        # Get miller indices from saved reflectons
        miller_indices = [r.miller_index for r in self._saved_reflections]

        self.miller_indices = flex.miller_index(miller_indices)

        print "Predicting new reflections"
        self.predict_reflections()

        # Put coords from same hkl in dict for saved reflections
        coord1 = defaultdict(list)
        for r1 in self._saved_reflections:
            coord1[r1.miller_index].append(r1.image_coord_px)

        # Put coords from same hkl in dict for new reflections
        coord2 = defaultdict(list)
        for r2 in self._new_reflections:
            coord2[r2.miller_index].append(r2.image_coord_px)

        # Print out coords for each hkl
        for h in coord1.keys():
            c1 = coord1[h]
            c2 = coord2[h]
            #print c1, c2

    def predict_reflections(self):
        '''Predict all reflection positions after refinement and make the
        bounding boxes.'''

        from dials.algorithms.spot_prediction import IndexGenerator

        from dials.algorithms.spot_prediction import ray_intersection
        from dials.algorithms.spot_prediction import reflection_frames
        from dials.algorithms.shoebox import BBoxCalculator
        from math import pi
        from dials.algorithms.refinement.prediction.predictors import \
                ScanVaryingReflectionListGenerator

        s0 = self.beam.get_s0()
        dmin = self.detector.get_max_resolution_at_corners(s0,
                                            self.beam.get_wavelength())
        sv_predictor = ScanVaryingReflectionListGenerator(self.pred_param,
                            self.beam, self.gonio, self.scan, dmin)

        # Duck typing to determine whether prediction is scan-varying or not
        try:
            refs = sv_predictor()
            self._new_reflections = ray_intersection(self.detector, refs)

        except AttributeError: # prediction seems to be scan-static
            from dials.algorithms.spot_prediction import RayPredictor
            m2 = self.gonio.get_rotation_axis()
            UB = self.crystal.get_U() * self.crystal.get_B()
            dphi = self.scan.get_oscillation_range(deg=False)
            predict_rays = RayPredictor(s0, m2, dphi)
            miller_indices = flex.miller_index(
                [r.miller_index for r in self._saved_reflections])
            self._new_reflections = reflection_frames(self.scan, ray_intersection(
                self.detector, predict_rays(miller_indices, UB)))

        self.sigma_divergence = self.beam.get_sigma_divergence()
        self.sigma_mosaicity = self.crystal.get_mosaicity()

        # Set the divergence and mosaicity
        n_sigma = 5.0
        delta_divergence = n_sigma * self.sigma_divergence * pi / 180.0
        delta_mosaicity = n_sigma * self.sigma_mosaicity * pi / 180.0

        # FIXME: DIALS_ASSERT(delta_divergence > 0.0) failure.
        #
        # Create the bounding box calculator
        #calculate_bbox = BBoxCalculator(self.beam, self.detector, self.gonio,
        #    self.scan, delta_divergence, delta_mosaicity)

        # Calculate the frame numbers of all the reflections
        #calculate_bbox(self._new_reflections)

        return self._new_reflections


class RefinerFactory(object):
    ''' Factory class to create refiners '''

    @staticmethod
    def from_parameters(params, verbosity):
        ''' Given a set of parameters, construct the refiner

        Params:
            params The input parameters
            verbosity The verbosity level

        Returns:
            The refiner instance

        '''

        parameterisation_strategy = \
                    RefinerFactory.configure_parameterisation(params)
        refinery_strategy = RefinerFactory.configure_refinery(params)
        reflections_strategy = RefinerFactory.configure_refman(params)
        target_strategy = RefinerFactory.configure_target(params)

        return Refiner(parameterisation_strategy, refinery_strategy,
                 reflections_strategy, target_strategy, verbosity)

    @staticmethod
    def configure_parameterisation(params):
        ''' Given a set of parameters, configure a factory to build a
        parameterisation from a set of experimental models

        Params:
            params The input parameters

        Returns:
            The parameterisation factory instance
        '''

        # Shorten parameter paths
        beam_options = params.refinement.parameterisation.beam
        crystal_options = params.refinement.parameterisation.crystal
        detector_options = params.refinement.parameterisation.detector

        return ParameterisationFactory(beam_options, crystal_options,
                                       detector_options)

    @staticmethod
    def configure_refinery(params):
        ''' Given a set of parameters, configure a factory to build a
        refinery

        Params:
            params The input parameters

        Returns:
            The refinery factory instance
        '''

        # Shorten parameter path
        options = params.refinement.refinery
        return RefineryFactory(options)

    @staticmethod
    def configure_refman(params):
        ''' Given a set of parameters, configure a factory to build a
        reflection manager

        Params:
            params The input parameters

        Returns:
            The reflection manager factory instance
        '''

        # Shorten parameter path
        options = params.refinement.reflections
        return RefmanFactory(options)

    @staticmethod
    def configure_target(params):
        ''' Given a set of parameters, configure a factory to build a
        target function

        Params:
            params The input parameters

        Returns:
            The target factory instance
        '''

        # Shorten parameter path
        options = params.refinement.target
        return TargetFactory(options)


class ParameterisationFactory(object):

    def __init__(self, beam_options, crystal_options, detector_options):

        # Shorten paths
        import dials.algorithms.refinement.parameterisation as par

        # Beam
        self._beam_par = par.BeamParameterisationOrientation
        self._beam_fix_all = beam_options.fix_beam
        self._beam_fix_in_spindle_plane = beam_options.fix_in_spindle_plane

        # Crystal
        self._crystal_fix_cell = crystal_options.fix_cell
        self._crystal_fix_orientation = crystal_options.fix_orientation
        self._crystal_scan_varying = crystal_options.scan_varying
        self._crystal_num_intervals = crystal_options.num_intervals

        if self._crystal_scan_varying:
            cop = par.ScanVaryingCrystalOrientationParameterisation
            cucp = par.ScanVaryingCrystalUnitCellParameterisation
        else:
            cop = par.CrystalOrientationParameterisation
            cucp = par.CrystalUnitCellParameterisation
        self._crystal_ori_par = cop
        self._crystal_uc_par = cucp

        # Detector
        if detector_options.panels == "single":
            dp = par.DetectorParameterisationSinglePanel
        elif detector_options.panels == "multiple":
            raise RuntimeError, "multiple panel detector parameterisation is not yet implemented"
        else:
            raise RuntimeError, "detector parameterisation type not recognised"

        self.detector_par = dp
        self._fix_detector = detector_options.fix_detector

        # Prediction equation parameterisation
        if self._crystal_scan_varying:
            pep = par.VaryingCrystalPredictionParameterisation
        else:
            pep = par.DetectorSpacePredictionParameterisation

        self.prediction_par = pep

    def __call__(self, beam, crystal, goniometer, detector, scan):

        beam_param = self._beam_par(beam, goniometer)
        if self._beam_fix_in_spindle_plane:
            beam_param.set_fixed([True, False])
        if self._beam_fix_all:
            beam_param.set_fixed([True, True])

        if self._crystal_scan_varying:
            xl_ori_param = self._crystal_ori_par(crystal,
                                                 scan.get_image_range(),
                                                 self._crystal_num_intervals)
            xl_uc_param = self._crystal_uc_par(crystal,
                                               scan.get_image_range(),
                                               self._crystal_num_intervals)
        else:
            xl_ori_param = self._crystal_ori_par(crystal)
            xl_uc_param = self._crystal_uc_par(crystal)

        if self._crystal_fix_orientation:
            xl_ori_param.set_fixed([True] * xl_ori_param.num_free())
        if self._crystal_fix_cell:
            xl_uc_param.set_fixed([True] * xl_uc_param.num_free())

        det_param = self.detector_par(detector)
        if self._fix_detector:
            det_param.set_fixed([True] * det_param.num_free())

        pred_param = self.prediction_par(detector, beam, crystal, goniometer,
                [det_param], [beam_param], [xl_ori_param], [xl_uc_param])

        return (beam_param, xl_ori_param, xl_uc_param, det_param, pred_param)

class RefineryFactory(object):

    def __init__(self, options):

        import dials.algorithms.refinement.engine as engine

        if options.engine  == "SimpleLBFGS":
            from engine import SimpleLBFGS as ref
        elif options.engine == "LBFGScurvs":
            from engine import LBFGScurvs as ref
        elif options.engine == "GaussNewtonIterations":
            from engine import GaussNewtonIterations as ref
        else:
            raise RuntimeError, "Refinement engine " + options.engine + " not recognised"

        self._refinery = ref
        self._track_step = options.track_step
        self._track_gradient = options.track_gradient
        self._logfile = options.log

    def __call__(self, target, prediction_parameterisation, verbosity):

        return self._refinery(
            target = target,
            prediction_parameterisation = prediction_parameterisation,
            log = self._logfile,
            verbosity = verbosity,
            track_step = self._track_step,
            track_gradient = self._track_gradient)

class RefmanFactory(object):

    def __init__(self, options):

        import dials.algorithms.refinement.target as target

        self._ref_per_degree = options.reflections_per_degree
        if options.use_all_reflections:
            self._ref_per_degree = None

        self._min_num_obs = options.minimum_number_of_reflections

        from target import ReflectionManager as refman
        self._refman = refman

    def __call__(self, reflections, beam, goniometer, scan, verbosity):

        from scitbx import matrix
        from math import sqrt

        # pull out data needed for refinement
        temp = [(ref.miller_index, ref.entering, ref.frame_number,
                 ref.rotation_angle, matrix.col(ref.beam_vector),
                 ref.image_coord_mm, ref.centroid_variance) \
                    for ref in reflections]
        (hkls, enterings, frames, angles,
         svecs, intersects, variances) = zip(*temp)

        # tease apart tuples to separate lists
        d1s, d2s = zip(*intersects)
        var_d1s, var_d2s, var_angles = zip(*variances)

        # change variances to sigmas
        sig_d1s = [sqrt(e) for e in var_d1s]
        sig_d2s = [sqrt(e) for e in var_d2s]
        sig_angles = [sqrt(e) for e in var_angles]

        assert len(hkls) == len(svecs) == len(d1s) == len(d2s) == \
               len(sig_d2s) == len(angles) == len(sig_angles)

        return self._refman(h_obs=hkls,
                            entering_obs=enterings,
                            frame_obs=frames,
                            svec_obs=svecs,
                            x_obs=d1s, sigx_obs=sig_d1s,
                            y_obs=d2s, sigy_obs=sig_d2s,
                            phi_obs=angles, sigphi_obs=sig_angles,
                            beam=beam,
                            gonio=goniometer,
                            scan=scan,
                            verbosity=verbosity,
                            nref_per_degree=self._ref_per_degree,
                            min_num_obs=self._min_num_obs)

class TargetFactory(object):

    def __init__(self, options):

        import dials.algorithms.refinement.target as target

        if options.implementation == "basic":
            from target import \
                LeastSquaresPositionalResidualWithRmsdCutoff as targ
        else:
            raise RuntimeError, "Target type " + options.implementation + \
                                " not recognised"

        self._frac_binsize_cutoff = options.bin_size_fraction
        if options.rmsd_cutoff == "fraction_of_bin_size":
            self._absolute_cutoffs = None
        elif options.rmsd_cutoff == "absolute":
            self._absolute_cutoffs = options.absolute_cutoffs
        else:
            raise RuntimeError, "Target function rmsd_cutoff option" + \
                options.rmsd_cutoff + " not recognised"

        # Reflection prediction
        from dials.algorithms.refinement.prediction import \
            ReflectionPredictor as rp

        self._target = targ
        self._ref_predictor = rp

    def __call__(self, crystal, beam, goniometer, detector, scan,
        refman, pred_param):

        sweep_range = scan.get_oscillation_range(deg=False)
        image_width = scan.get_oscillation(deg=False)[1]

        rp = self._ref_predictor(crystal, beam, goniometer, sweep_range)
        return self._target(rp, detector, refman, pred_param,
                            image_width, self._frac_binsize_cutoff,
                            self._absolute_cutoffs)
