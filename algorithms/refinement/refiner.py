#!/usr/bin/env python
#
# dials.algorithms.refinement.refiner.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class Refiner(object):
    ''' The refiner class. '''

    def __init__(self, engine, scan_varying, verbosity):
        ''' Initialise the refiner class. '''

        # Set the parameters
        self.engine = engine
        self.scan_varying = scan_varying
        self.verbosity = verbosity

    def __call__(self, sweep, crystal, reflections):
        ''' Call to refine.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflection list

        Returns:
            The refined (sweep, crystal)

        '''
        from scitbx import matrix
        from math import sqrt

        # Get the models from the sweep
        self.beam = sweep.get_beam()
        self.detector = sweep.get_detector()
        self.gonio = sweep.get_goniometer()
        self.scan = sweep.get_scan()
        self.crystal = crystal

        self.UB = self.crystal.get_U() * self.crystal.get_B()
        self.sigma_divergence = self.beam.get_sigma_divergence()
        self.sigma_mosaicity = self.crystal.get_mosaicity()

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

        # pull out data needed for refinement
        temp = [(ref.miller_index, ref.entering, ref.frame_number,
                 ref.rotation_angle, matrix.col(ref.beam_vector),
                 ref.image_coord_mm, ref.centroid_variance) \
                    for ref in self.reflections]
        hkls, enterings, frames, angles, svecs, intersects, variances = zip(*temp)

        # tease out tuples to separate lists
        d1s, d2s = zip(*intersects)
        var_d1s, var_d2s, var_angles = zip(*variances)

        px_size = self.detector.get_pixel_size()
        im_width = self.scan.get_oscillation(deg=False)[1]

        # get the angular range of the sweep
        sweep_range = self.scan.get_oscillation_range(deg=False)

        # change variances to sigmas
        sig_d1s = [sqrt(e) for e in var_d1s]
        sig_d2s = [sqrt(e) for e in var_d2s]
        sig_angles = [sqrt(e) for e in var_angles]

        # DEBUGGING: ignore calculated variances and just use invented values,
        # based on half the pixel size and half the image width
        #sig_d1s = [px_size[0] / 2.] * len(hkls)
        #sig_d2s = [px_size[1] / 2.] * len(hkls)
        #sig_angles = [im_width / 2.] * len(hkls)

        assert len(hkls) == len(svecs) == len(d1s) == len(d2s) == \
               len(sig_d2s) == len(angles) == len(sig_angles)

        from dials.algorithms.refinement import print_model_geometry, refine

        print "Prior to refinement the experimental model is:"
        print_model_geometry(self.beam, self.detector, self.crystal)

        refine(self.beam, self.gonio, self.crystal, self.detector, im_width,
               self.scan, hkls, enterings, frames,
               svecs, d1s, sig_d1s, d2s, sig_d2s, angles,
               sig_angles, verbosity = self.verbosity, fix_cell=False,
               scan_varying=self.scan_varying)

        print
        print "Refinement has completed with the following geometry:"
        print_model_geometry(self.beam, self.detector, self.crystal)

        # Do a test of new reflection pos
        self._update_reflections_test()

        # Update models in sweep
        sweep.set_detector(self.detector)
        sweep.set_beam(self.beam)
        sweep.set_goniometer(self.gonio)

        # Return refined models
        return sweep, crystal

    def _update_reflections_test(self):
        from cctbx.array_family import flex
        from collections import defaultdict

        # Get miller indices from saved reflectons
        miller_indices = [r.miller_index for r in self._saved_reflections]

        self.miller_indices = flex.miller_index(miller_indices)

        print "Predicting new reflections"
        self._predict_reflections()

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

    def _predict_reflections(self):
        '''Predict the reflection positions and bounding boxes.'''
        from dials.algorithms.spot_prediction import IndexGenerator
        from dials.algorithms.spot_prediction import RayPredictor
        from dials.algorithms.spot_prediction import ray_intersection
        from dials.algorithms.spot_prediction import reflection_frames
        from dials.algorithms.integration import BBoxCalculator
        from math import pi

        # Create the spot predictor
        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        UB = self.UB
        dphi = self.scan.get_oscillation_range(deg=False)
        predict_rays = RayPredictor(s0, m2, dphi)

        # Predict the reflections
        self._new_reflections = reflection_frames(self.scan, ray_intersection(
            self.detector, predict_rays(self.miller_indices, UB)))

        # Set the divergence and mosaicity
        n_sigma = 5.0
        delta_divergence = n_sigma * self.sigma_divergence * pi / 180.0
        delta_mosaicity = n_sigma * self.sigma_mosaicity * pi / 180.0

        # Create the bounding box calculator
        calculate_bbox = BBoxCalculator(self.beam, self.detector, self.gonio,
            self.scan, delta_divergence, delta_mosaicity)

        # Calculate the frame numbers of all the reflections
        calculate_bbox(self._new_reflections)

class RefinerFactory(object):
    ''' Factory class to create refiners '''

    @staticmethod
    def from_parameters(params, verbosity):
        ''' Given a set of parametets, construct the refiner

        Params:
            params The input parameters
            verbosity The verbosity level

        Returns:
            The refiner instance

        '''
        return Refiner(
            engine=params.refinement.engine,
            scan_varying=params.refinement.scan_varying,
            verbosity=verbosity)
