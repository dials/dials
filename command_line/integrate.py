#!/usr/bin/env python
#
# dials.integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.util.script import ScriptRunner


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


class Integrator(object):
    ''' The integrator base class. '''

    def __init__(self, compute_spots, compute_background, compute_intensity):
        ''' Initialise the integrator base class.

        Params:
            compute_spots The spot extractor strategy
            compute_background The background strategy
            compute_intensity The intensity strategy

        '''
        # Initialise the reflection extractor
        self.compute_spots = compute_spots
        self.compute_background = compute_background
        self.compute_intensity = compute_intensity

    def __call__(self, sweep, crystal, reflections=None):
        ''' Call to integrate.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflection list

        Returns:
            A reflection list

        '''
        # Extract the reflections from the sweep
        if reflections == None:
            reflections = self.compute_spots(sweep, crystal)

        # Calculate the background
        reflections = self.compute_background(sweep, crystal, reflections)

        # Calculate the intensity and return
        return self.compute_intensity(sweep, crystal, reflections)


class IntegratorFactory(object):
    ''' Factory class to create integrators '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parametets, construct the integrator

        Params:
            params The input parameters

        Returns:
            The integrator instance

        '''

        # Configure the algorithms to extract reflections, compute the
        # background intensity and integrate the reflection intensity
        compute_spots = IntegratorFactory.configure_extractor(params)
        compute_background = IntegratorFactory.configure_background(params)
        compute_intensity = IntegratorFactory.configure_intensity(params)

        # Return the integrator with the given strategies
        return Integrator(compute_spots=compute_spots,
                          compute_background=compute_background,
                          compute_intensity=compute_intensity)

    @staticmethod
    def configure_extractor(params):
        ''' Given a set of parameters, configure the reflection extractor

        Params:
            params The input parameters

        Returns:
            The extractor instance

        '''
        return ReflectionExtractor(params.integration.shoebox.n_sigma)

    @staticmethod
    def configure_background(params):
        ''' Given a set of parameters, configure the background calculator

        Params:
            params The input parameters

        Returns:
            The background calculator instance

        '''
        from dials.algorithms.background import NullSubtractor
        from dials.algorithms.background import XdsSubtractor
        from dials.algorithms.background import FableSubtractor

        # Configure the NULL subtractor
        if params.integration.background.algorithm == 'none':
            algorithm = NullSubtractor()

        # Configure the XDS subtractor
        elif params.integration.background.algorithm == 'xds':
            algorithm = XdsSubtractor(
                min_data=params.integration.background.min_pixels,
                n_sigma=params.integration.background.n_sigma)

        # Configure the Fable subtractor
        elif params.integration.background.algorithm == 'fable':
            algorithm = FableSubtractor(
                min_data=params.integration.background.min_pixels,
                n_sigma=params.integration.background.n_sigma)

        # Configure the flat subtractor
        elif params.integration.background.algorithm == 'flat':
            raise RuntimeError('Not implemented yet')

        # Configure the curved subtractor
        elif params.integration.background.algorithm == 'curved':
            raise RuntimeError('Not implemented yet')

        # Configure the esmerelda subtractor
        elif params.integration.background.algorithm == 'esmerelda':
            raise RuntimeError('Not implemented yet')

        # Unknown subtractor
        else:
            raise RuntimeError('Unknown background algorithm')

        # Return the algorithm
        return algorithm

    @staticmethod
    def configure_intensity(params):
        ''' Given a set of parameters, configure the intensity calculator

        Params:
            params The input parameters

        Returns:
            The intensity calculator instance

        '''
        return lambda x, y, z: z


class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # The script usage
        usage = "usage: %prog [options] [param.phil] sweep.json crystal.json"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'integrated.pickle',
            help = 'Set the filename for integrated reflections.')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.model.serialize import load
        import cPickle as pickle

        # Check the number of arguments is correct
        if len(args) != 2:
            self.config().print_help()
            return

        # Get the integrator from the input parameters
        print 'Configurating integrator from input parameters'
        integrate = IntegratorFactory.from_parameters(params)

        # Try to load the models
        print 'Loading models from {0} and {1}'.format(args[0], args[1])
        sweep = load.sweep(args[0])
        crystal = load.crystal(args[1])

        # Intregate the sweep's reflections
        print 'Integrating reflections'
        reflections = integrate(sweep, crystal)

        # Save the reflections to file
        print 'Saving reflections to {0}'.format(options.output_filename)
        pickle.dump(reflections, open(options.output_filename, 'w'))


if __name__ == '__main__':
    script = Script()
    script.run()
