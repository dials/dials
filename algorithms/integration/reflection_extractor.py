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

    def __init__(self, bbox_nsigma, filter_by_bbox=False, filter_by_zeta=0.0,
                 filter_by_xds_small_angle=False, filter_by_xds_angle=False,
                 gain_map=None, dark_map=None, kernel_size=None,
                 n_sigma_b=None, n_sigma_s=None):
        ''' Initialise the extractor

        Params:
            bbox_nsigma The number of standard deviations for bbox
            filter_by_bbox True/False
            filter_by_zeta Minimum zeta value
            filter_by_xds_small_angle True/False
            filter_by_xds_angle True/False
            gain_map The detector gain map
            dark_map The detector dark map
            kernel_size The size of the smoothing kernel to use
            n_sigma_b The number of sigmas for the background
            n_sigma_s The number of sigmas for a strong pixel

        '''
        # Get parameters we need
        self.bbox_nsigma = bbox_nsigma
        self.filter_by_bbox = filter_by_bbox
        self.filter_by_zeta = filter_by_zeta
        self.filter_by_xds_small_angle = filter_by_xds_small_angle
        self.filter_by_xds_angle = filter_by_xds_angle
        self.gain_map = gain_map
        self.dark_map = dark_map
        self.kernel_size = kernel_size
        self.n_sigma_b = n_sigma_b
        self.n_sigma_s = n_sigma_s

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
        from dials.algorithms.shoebox import BBoxCalculator

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
        from dials.algorithms import shoebox
        from dials.algorithms.integration import extract_reflection_profiles
        from dials.algorithms.integration import filter_by_detector_mask
        from dials.algorithms.integration import filter
        from math import sqrt

        # Get models from the sweep
        beam = sweep.get_beam()
        detector = sweep.get_detector()
        gonio = sweep.get_goniometer()
        scan = sweep.get_scan()

        # Get the detector mask
        detector_mask = sweep[0] >= 0

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
        overlaps = shoebox.find_overlapping(reflections)
        Command.end('Found {0} overlaps'.format(len(overlaps)))

        # Set all reflections which overlap bad pixels to zero
        Command.start('Filtering reflections by detector mask')
        array_range = scan.get_array_range()
        filter_by_detector_mask(reflections, detector_mask, array_range)
        Command.end('Filtered {0} reflections by detector mask'.format(
            len([r for r in reflections if r.is_valid()])))

        # Filter the reflections by bbox
        if self.filter_by_bbox:
            Command.start('Filtering reflections by bounding box')
            filter.by_bbox_volume(reflections)
            Command.end('Filtered {0} reflections by bounding box'.format(
                len([r for r in reflections if r.is_valid()])))

        # Filter the reflections by zeta
        if self.filter_by_zeta > 0:
            Command.start('Filtering reflections by zeta >= {0}'.format(
                self.filter_by_zeta))
            filter.by_zeta(gonio, beam, reflections, self.filter_by_zeta)
            Command.end('Filtered {0} reflections by zeta >= {1}'.format(
                len([r for r in reflections if r.is_valid()]),
                self.filter_by_zeta))

        # Filter the reflections by xds small angle validity
        if self.filter_by_xds_small_angle:
            Command.start('Filtering reflections by xds small angle')
            filter.by_xds_small_angle(gonio, beam, reflections,
              self.bbox_nsigma * crystal.get_mosaicity(deg=False))
            Command.end('Filtered {0} reflections by xds small angle'.format(
                len([r for r in reflections if r.is_valid()])))

        # Filter the reflections by xds small angle validity
        if self.filter_by_xds_angle:
            Command.start('Filtering reflections by xds angle')
            filter.by_xds_angle(gonio, beam, reflections,
              self.bbox_nsigma * crystal.get_mosaicity(deg=False))
            Command.end('Filtered {0} reflections by xds angle'.format(
                len([r for r in reflections if r.is_valid()])))

        # Extract the reflection profiles
        extract_reflection_profiles(
            sweep, reflections, overlaps,
            self.gain_map, self.dark_map,
            self.kernel_size,
            self.n_sigma_b, self.n_sigma_s,
            detector_mask)

        # Return the list of reflections
        return reflections
