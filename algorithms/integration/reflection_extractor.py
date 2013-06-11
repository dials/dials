#!/usr/bin/env python
#
# dials.algorithms.integration.reflection_extractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class ReflectionExtractor(object):
    ''' Class to extract basic reflection information. '''

    def __init__(self, bbox_nsigma):
        ''' Initialise the extractor

        Params:
            bbox_nsigma The number of standard deviations for bbox

        '''
        # Get parameters we need
        self.bbox_nsigma = bbox_nsigma

    def __call__(self, sweep, crystal):
        ''' Extract the basic reflection properties from the sweep

        Params:
            sweep The sweep to process
            crystal The crystal model to use

        Returns:
            A list of reflections

        '''
        # Initialise the algorithms
        self.initialise(sweep, crystal)

        # Extract the reflections
        return self.extract(sweep, crystal)

    def initialise(self, sweep, crystal):
        ''' Initialise the extraction algorithms

        Params:
            sweep The sweep to process
            crystal The crystal model to use

        '''
        from cctbx import sgtbx
        from dials.util.command_line import Command
        from dials.algorithms.spot_prediction import IndexGenerator
        from dials.algorithms.spot_prediction import RayPredictor
        from dials.algorithms.integration import BBoxCalculator

        # Get models from the sweep
        beam = sweep.get_beam()
        detector = sweep.get_detector()
        gonio = sweep.get_goniometer()
        scan = sweep.get_scan()

        # Create the index generator
        self.generate_hkl = IndexGenerator(
            crystal.get_unit_cell(),
            sgtbx.space_group_type(crystal.get_space_group()),
            detector.get_max_resolution(beam.get_s0(), beam.get_wavelength()))

        # Create the spot predictor
        self.predict_rays = RayPredictor(
            beam.get_s0(),
            gonio.get_rotation_axis(),
            scan.get_oscillation_range(deg=False))

        # Create the bbox calculator
        self.compute_bbox = BBoxCalculator(
            beam, detector, gonio, scan,
            self.bbox_nsigma * beam.get_sigma_divergence(deg=False),
            self.bbox_nsigma * crystal.get_mosaicity(deg=False))

    def extract(self, sweep, crystal):
        ''' Extract the reflections from the sweep.

        Params:
            sweep The sweep to process
            crystal The crystal model to use

        Returns:
            A list of reflections

        '''
        from dials.util.command_line import Command
        from dials.algorithms.spot_prediction import ray_intersection
        from dials.algorithms.spot_prediction import reflection_frames
        from dials.algorithms.integration import find_overlapping_reflections
        from dials.algorithms.integration import extract_reflection_profiles

        # Get models from the sweep
        beam = sweep.get_beam()
        detector = sweep.get_detector()
        gonio = sweep.get_goniometer()
        scan = sweep.get_scan()

        # Generate Indices
        Command.start('Generating miller indices')
        miller_indices = self.generate_hkl.to_array()
        Command.end('Generating {0} miller indices'.format(len(miller_indices)))

        # Predict reflections
        Command.start('Predicting rays')
        reflections = self.predict_rays(miller_indices, crystal.get_A())
        Command.end('Predicted {0} rays'.format(len(reflections)))

        # Get detector coordinates (mm)
        Command.start('Calculating ray-detector intersections')
        reflections = ray_intersection(detector, reflections)
        Command.end('Calculated {0} intersections'.format(len(reflections)))

        # Calculate the frame numbers of all the reflections
        Command.start('Calculating reflection frames')
        reflections = reflection_frames(scan, reflections)
        Command.end('Calculated {0} frames'.format(len(reflections)))

        # Calculate the bounding boxes of all the reflections
        Command.start('Calculating bounding boxes')
        self.compute_bbox(reflections)
        Command.end('Calculated {0} bounding boxes'.format(len(reflections)))

        # Find overlapping reflections
        Command.start('Finding overlapping reflections')
        overlaps = find_overlapping_reflections(reflections)
        Command.end('Found {0} overlaps'.format(len(overlaps)))

        # Extract the reflection profiles
        extract_reflection_profiles(sweep, reflections, overlaps)

        # Return the list of reflections
        return reflections
