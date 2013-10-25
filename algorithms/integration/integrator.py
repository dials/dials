#!/usr/bin/env python
#
# dials.algorithms.integration.integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class Integrator(object):
    ''' The integrator base class. '''

    def __init__(self, compute_spots, compute_background, compute_centroid,
                       compute_intensity, correct_intensity):
        ''' Initialise the integrator base class.

        Params:
            compute_spots The spot extractor strategy
            compute_background The background strategy
            compute_centroid The centroid strategy
            compute_intensity The intensity strategy
            correct_intensity The intensity correction strategy

        '''
        self.compute_spots = compute_spots
        self.compute_background = compute_background
        self.compute_centroid = compute_centroid
        self.compute_intensity = compute_intensity
        self.correct_intensity = correct_intensity

    def __call__(self, sweep, crystal, reflections = None, reference = None):
        ''' Call to integrate.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflection list
            reference The reference profiles

        Returns:
            A reflection list

        '''

        # Extract the reflections from the sweep
        reflections = self.compute_spots(sweep, crystal, reflections)

        # Process the reflections in the following way:
        #  - First calculate the background
        #  - Second calculate the reflection centroids
        #  - Third calculate the reflection intensities
        #  - Fourth correct the reflection intensities
        reflections = self.compute_background(sweep, crystal, reflections)
        reflections = self.compute_centroid(sweep, crystal, reflections)
        reflections = self.compute_intensity(sweep, crystal, reflections, reference)
        reflections = self.correct_intensity(sweep, crystal, reflections)

        # Return the reflections
        return reflections


class Integrator2(object):
    ''' The integrator base class. '''

    def __init__(self, n_sigma, n_blocks, filter_by_zeta,
                 compute_background, compute_centroid,
                 compute_intensity, correct_intensity):
        ''' Initialise the integrator base class.

        Params:
            compute_spots The spot extractor strategy
            compute_background The background strategy
            compute_centroid The centroid strategy
            compute_intensity The intensity strategy
            correct_intensity The intensity correction strategy

        '''
        self.n_sigma = n_sigma
        self.n_blocks = n_blocks
        self.filter_by_zeta = filter_by_zeta
        self.compute_background = compute_background
        self.compute_centroid = compute_centroid
        self.compute_intensity = compute_intensity
        self.correct_intensity = correct_intensity

    def __call__(self, sweep, crystal, reference=None, extracted=None):
        ''' Call to integrate.

        Params:
            sweep The sweep to process
            crystal The crystal to process
            reflections The reflection list
            reference The reference profiles

        Returns:
            A reflection list

        '''
        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.integration import ReflectionBlockExtractor
        from dials.model.data import ReflectionList
        from dials.array_family import flex
        from dials.util.command_line import Command

        # Predict a load of reflections
        if extracted == None:
            predict = ReflectionPredictor()
            predicted = predict(sweep, crystal)
        else:
            predicted = None

        # Get the extractor
        extract = ReflectionBlockExtractor(sweep, crystal, predicted,
            self.n_sigma, self.n_blocks, self.filter_by_zeta, extracted)

        # Loop through all the blocks
        result = ReflectionList()
        print ''
        for reflections in extract:

            self.compute_background(sweep, crystal, reflections)
            self.compute_centroid(sweep, crystal, reflections)
            self.compute_intensity(sweep, crystal, reflections, reference)
            self.correct_intensity(sweep, crystal, reflections)

            # Remove the profiles from the reflections
            for r in reflections:
                r.shoebox = flex.double(flex.grid(0, 0, 0))
                r.shoebox_mask = flex.int(flex.grid(0, 0, 0))
                r.shoebox_background = flex.double(flex.grid(0, 0, 0))
                r.transformed_shoebox = flex.double(flex.grid(0, 0, 0))
                r.transformed_shoebox_background = flex.double(flex.grid(0, 0, 0))
            result.extend(reflections)
            print ''

        # Return the reflections
        return ReflectionList(sorted(result, key=lambda x: x.miller_index))


class IntegratorFactory(object):
    ''' Factory class to create integrators '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the integrator

        Params:
            params The input parameters

        Returns:
            The integrator instance

        '''
        # Configure the algorithms to extract reflections, compute the
        # background intensity and integrate the reflection intensity
        compute_spots = IntegratorFactory.configure_extractor(params)
        compute_background = IntegratorFactory.configure_background(params)
        compute_centroid = IntegratorFactory.configure_centroid(params)
        compute_intensity = IntegratorFactory.configure_intensity(params)
        correct_intensity = IntegratorFactory.configure_correction(params)

        # Return the integrator with the given strategies
        return Integrator(compute_spots = compute_spots,
                          compute_background = compute_background,
                          compute_centroid = compute_centroid,
                          compute_intensity = compute_intensity,
                          correct_intensity = correct_intensity)

    @staticmethod
    def configure_extractor(params):
        ''' Given a set of parameters, configure the reflection extractor

        Params:
            params The input parameters

        Returns:
            The extractor instance

        '''
        from dials.algorithms.integration import ReflectionExtractor

        # Load some lookup maps
        mask = IntegratorFactory.load_image(params.lookup.mask)
        gain = IntegratorFactory.load_image(params.lookup.gain_map)
        dark = IntegratorFactory.load_image(params.lookup.dark_map)

        # Shorten parameter path
        integration = params.integration

        # Return the reflection extractor instance
        return ReflectionExtractor(
            bbox_nsigma = integration.shoebox.n_sigma,
            filter_by_bbox = integration.filter.by_bbox,
            filter_by_zeta = integration.filter.by_zeta,
            filter_by_xds_small_angle = integration.filter.by_xds_small_angle,
            filter_by_xds_angle = integration.filter.by_xds_angle,
            mask = mask,
            gain = gain,
            dark = dark,
            kernel_size = integration.shoebox.kernel_size,
            n_sigma_b = integration.shoebox.sigma_background,
            n_sigma_s = integration.shoebox.sigma_strong)

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
        from dials.algorithms.background import FlatSubtractor
        from dials.algorithms.background import CurvedSubtractor
        from dials.algorithms.background import InclinedSubtractor

        # Shorten parameter path
        integration = params.integration

        # Configure the NULL subtractor
        if (integration.background.algorithm == 'none' or
            integration.background.algorithm == None):
            algorithm = NullSubtractor()

        # Configure the XDS subtractor
        elif integration.background.algorithm == 'xds':
            algorithm = XdsSubtractor(
                min_data = integration.background.min_pixels,
                n_sigma = integration.background.n_sigma)

        # Configure the Fable subtractor
        elif params.integration.background.algorithm == 'fable':
            algorithm = FableSubtractor(
                min_data = integration.background.min_pixels,
                n_sigma = integration.background.n_sigma)

        # Configure the flat subtractor
        elif integration.background.algorithm == 'flat':
            algorithm = FlatSubtractor()

        # Configure the inclined plane subtractor
        elif integration.background.algorithm == 'inclined':
            algorithm = InclinedSubtractor()

        # Configure the esmerelda curved subtractor
        elif integration.background.algorithm == 'esmeralda':
            algorithm = CurvedSubtractor()

        # Unknown subtractor
        else:
            raise RuntimeError('Unknown background algorithm')

        # Return the algorithm
        return algorithm

    @staticmethod
    def configure_centroid(params):
        ''' Given a set of parameters, configure the centroid calculator

        Params:
            params The input parameters

        Returns:
            The centroid calculator instance

        '''
        from dials.algorithms.centroid.centroider import Centroider
        return Centroider()

    @staticmethod
    def configure_intensity(params):
        ''' Given a set of parameters, configure the intensity calculator

        Params:
            params The input parameters

        Returns:
            The intensity calculator instance

        '''
        from dials.algorithms.integration import Summation2d
        from dials.algorithms.integration import Summation3d
        from dials.algorithms.integration import SummationReciprocalSpace
        from dials.algorithms.integration import ProfileFittingReciprocalSpace
        from dials.algorithms.integration.mosflm_like import MosflmProfileFitting

        # Shorten parameter path
        integration = params.integration

        # Configure the 2D summation algorithm
        if integration.algorithm == 'sum2d':
            algorithm = Summation2d()

        # Configure the 3D summation algorithm
        elif integration.algorithm == 'sum3d':
            algorithm = Summation3d()

        # Configure the reciprocal space summation algorithm
        elif params.integration.algorithm == 'sum_rs':
            algorithm = SummationReciprocalSpace(
                n_sigma = integration.shoebox.n_sigma,
                grid_size = integration.reciprocal_space.grid_size)

        # Configure the 2D profile fitting algorithm
        elif integration.algorithm == 'fit_2d':
            algorithm = MosflmProfileFitting(
                nblocks = integration.mosflm.nblocks)

        # Configure the 3D profile fitting algorithm
        elif integration.algorithm == 'fit_3d':
            raise RuntimeError('Not implemented yet')

        # Configure the reciprocal space profile fitting algorithm
        elif integration.algorithm == 'fit_rs':
            algorithm = ProfileFittingReciprocalSpace(
                n_sigma = integration.shoebox.n_sigma,
                grid_size = integration.reciprocal_space.grid_size,
                frame_interval = integration.profile.reference_frame_interval,
                threshold = integration.profile.reference_signal_threshold)

        # Unknown algorithm
        else:
            raise RuntimeError('Unknown integration algorithm')

        # Return the algorithm
        return algorithm

    @staticmethod
    def configure_correction(params):
        ''' Given a set of parameters, configure the intensity correction

        Params:
            params The input parameters

        Returns:
            The intensity corrector instance

        '''
        from dials.algorithms.integration.lp_correction import correct_intensity
        return correct_intensity

    @staticmethod
    def load_image(filename):
        ''' Given a filename, load an image

        Params:
            filename The input filename

        Returns:
            The image or None

        '''
        from dials.util import image

        # If no filename is set then return None
        if not filename:
            return None

        # Read the image and return the image data
        handle = image.reader()
        handle.read_file(filename)
        return handle.get_data()

class IntegratorFactory2(object):
    ''' Factory class to create integrators '''

    @staticmethod
    def from_parameters(params):
        ''' Given a set of parameters, construct the integrator

        Params:
            params The input parameters

        Returns:
            The integrator instance

        '''
        # Configure the algorithms to extract reflections, compute the
        # background intensity and integrate the reflection intensity
        compute_background = IntegratorFactory2.configure_background(params)
        compute_centroid = IntegratorFactory2.configure_centroid(params)
        compute_intensity = IntegratorFactory2.configure_intensity(params)
        correct_intensity = IntegratorFactory2.configure_correction(params)

        # Return the integrator with the given strategies
        return Integrator2(n_sigma = params.integration.shoebox.n_sigma,
                          n_blocks = params.integration.shoebox.n_blocks,
                          filter_by_zeta = params.integration.filter.by_zeta,
                          compute_background = compute_background,
                          compute_centroid = compute_centroid,
                          compute_intensity = compute_intensity,
                          correct_intensity = correct_intensity)

    @staticmethod
    def configure_extractor(params):
        ''' Given a set of parameters, configure the reflection extractor

        Params:
            params The input parameters

        Returns:
            The extractor instance

        '''
        from dials.algorithms.integration import ReflectionExtractor

        # Load some lookup maps
        mask = IntegratorFactory.load_image(params.lookup.mask)
        gain = IntegratorFactory.load_image(params.lookup.gain_map)
        dark = IntegratorFactory.load_image(params.lookup.dark_map)

        # Shorten parameter path
        integration = params.integration

        # Return the reflection extractor instance
        return ReflectionExtractor(
            bbox_nsigma = integration.shoebox.n_sigma,
            filter_by_bbox = integration.filter.by_bbox,
            filter_by_zeta = integration.filter.by_zeta,
            filter_by_xds_small_angle = integration.filter.by_xds_small_angle,
            filter_by_xds_angle = integration.filter.by_xds_angle,
            mask = mask,
            gain = gain,
            dark = dark,
            kernel_size = integration.shoebox.kernel_size,
            n_sigma_b = integration.shoebox.sigma_background,
            n_sigma_s = integration.shoebox.sigma_strong)

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
        from dials.algorithms.background import FlatSubtractor
        from dials.algorithms.background import CurvedSubtractor
        from dials.algorithms.background import InclinedSubtractor

        # Shorten parameter path
        integration = params.integration

        # Configure the NULL subtractor
        if (integration.background.algorithm == 'none' or
            integration.background.algorithm == None):
            algorithm = NullSubtractor()

        # Configure the XDS subtractor
        elif integration.background.algorithm == 'xds':
            algorithm = XdsSubtractor(
                min_data = integration.background.min_pixels,
                n_sigma = integration.background.n_sigma)

        # Configure the Fable subtractor
        elif params.integration.background.algorithm == 'fable':
            algorithm = FableSubtractor(
                min_data = integration.background.min_pixels,
                n_sigma = integration.background.n_sigma)

        # Configure the flat subtractor
        elif integration.background.algorithm == 'flat':
            algorithm = FlatSubtractor()

        # Configure the inclined plane subtractor
        elif integration.background.algorithm == 'inclined':
            algorithm = InclinedSubtractor()

        # Configure the esmerelda curved subtractor
        elif integration.background.algorithm == 'esmeralda':
            algorithm = CurvedSubtractor()

        # Unknown subtractor
        else:
            raise RuntimeError('Unknown background algorithm')

        # Return the algorithm
        return algorithm

    @staticmethod
    def configure_centroid(params):
        ''' Given a set of parameters, configure the centroid calculator

        Params:
            params The input parameters

        Returns:
            The centroid calculator instance

        '''
        from dials.algorithms.centroid.centroider import Centroider
        return Centroider()

    @staticmethod
    def configure_intensity(params):
        ''' Given a set of parameters, configure the intensity calculator

        Params:
            params The input parameters

        Returns:
            The intensity calculator instance

        '''
        from dials.algorithms.integration import Summation2d
        from dials.algorithms.integration import Summation3d
        from dials.algorithms.integration import SummationReciprocalSpace
        from dials.algorithms.integration import ProfileFittingReciprocalSpace2
        from dials.algorithms.integration.mosflm_like import MosflmProfileFitting

        # Shorten parameter path
        integration = params.integration

        # Configure the 2D summation algorithm
        if integration.algorithm == 'sum2d':
            algorithm = Summation2d()

        # Configure the 3D summation algorithm
        elif integration.algorithm == 'sum3d':
            algorithm = Summation3d()

        # Configure the reciprocal space summation algorithm
        elif params.integration.algorithm == 'sum_rs':
            algorithm = SummationReciprocalSpace(
                n_sigma = integration.shoebox.n_sigma,
                grid_size = integration.reciprocal_space.grid_size)

        # Configure the 2D profile fitting algorithm
        elif integration.algorithm == 'fit_2d':
            algorithm = MosflmProfileFitting(
                nblocks = integration.mosflm.nblocks)

        # Configure the 3D profile fitting algorithm
        elif integration.algorithm == 'fit_3d':
            raise RuntimeError('Not implemented yet')

        # Configure the reciprocal space profile fitting algorithm
        elif integration.algorithm == 'fit_rs':

            from dials.algorithms.integration import ReferenceProfileFactory
            learner = ReferenceProfileFactory.from_parameters2(params)

            algorithm = ProfileFittingReciprocalSpace2(
                learner = learner,
                n_sigma = integration.shoebox.n_sigma,
                grid_size = integration.reciprocal_space.grid_size,
                frame_interval = integration.profile.reference_frame_interval,
                threshold = integration.profile.reference_signal_threshold)

        # Unknown algorithm
        else:
            raise RuntimeError('Unknown integration algorithm')

        # Return the algorithm
        return algorithm

    @staticmethod
    def configure_correction(params):
        ''' Given a set of parameters, configure the intensity correction

        Params:
            params The input parameters

        Returns:
            The intensity corrector instance

        '''
        from dials.algorithms.integration.lp_correction import correct_intensity
        return correct_intensity

    @staticmethod
    def load_image(filename):
        ''' Given a filename, load an image

        Params:
            filename The input filename

        Returns:
            The image or None

        '''
        from dials.util import image

        # If no filename is set then return None
        if not filename:
            return None

        # Read the image and return the image data
        handle = image.reader()
        handle.read_file(filename)
        return handle.get_data()
