#
# spot_finder.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.interfaces.peak_finding import SpotFinderInterface
from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy

class SpotFinder(SpotFinderInterface):
    '''A class to perform spot finding operations on a sweep of images.'''

    def __init__(self, **kwargs):
        '''Initialise the algorithm with some parameters.

        Params:
            kwargs The keyword arguments. This algorithm expects the following
            arguments:
                min_spot_size The minimum number of pixels in spot
                max_separation The maximum maximum-centroid distance

        '''
        # Get some parameters
        self._min_spot_size = kwargs.get('min_spot_size', 6)
        self._max_separation = kwargs.get('max_separation', 2)

        # Set the threshold strategy
        self._threshold_strategy = kwargs.get(
            'threshold_strategy', XDSThresholdStrategy())

    def __call__(self, sweep):
        '''The main function of the spot finder. Select the pixels from
        the sweep and then group the pixels into spots. Return the data
        in the form of a reflection list.

        Params:
            sweep The sweep object

        Returns:
            The reflection list

        '''
        from dials.util.command_line import Command

        # Set a command indent to 4
        Command.indent = 4
        print '\nFinding spots in {0} images...'.format(len(sweep))

        # Extract the image pixels from the sweep
        coords, intensity = self._extract_pixels(sweep)

        # Label the pixels and group into spots
        Command.start('Labelling connected components')
        labels = self._label_pixels(coords, sweep)
        Command.end('Found {0} connected components'.format(max(labels) + 1))

        # Filter spots that are too small
        Command.start('Filtering spots by size')
        spots = self._filter_spots(labels)
        Command.end('Filtered {0} spots by size'.format(len(spots)))

        # Calculate the bounding box for each spot
        Command.start('Calculating bounding boxes')
        bbox = self._calculate_bbox(coords, spots, sweep)
        Command.end('Calculated {0} bounding boxes'.format(len(bbox)))

        # Calculate the spot centroids
        Command.start('Calculating centroids')
        cpos, cvar, cerr, ctot = self._calculate_centroids(
            coords, intensity, spots)
        Command.end('Calculated {0} centroids'.format(len(cpos)))

        # Filter the spots by centroid-maxmimum distance
        Command.start('Filtering spots by distance')
        index = self._filter_maximum_centroid(coords, intensity, spots, cpos)
        Command.end('Filtered {0} spots by distance'.format(len(index)))

        # Create a reflection list and return
        return self._create_reflection_list(coords, intensity, spots, bbox,
            cpos, cvar, cerr, ctot, index)

    def _extract_pixels(self, sweep):
        '''Extract the pixels from the sweep

        Params:
            sweep The sweep object

        Returns:
            The list of selected pixels

        '''
        from dials.util.command_line import ProgressBar
        from scitbx.array_family import flex
        from dials.algorithms.peak_finding import flex_vec3_int

        # Initialise the pixel arrays
        coords = flex_vec3_int()
        intensity = flex.int()

        # Get the start index and trusted range from the sweep
        start = sweep.get_array_range()[0]
        trusted_range = sweep.get_detector().get_trusted_range()

        # Loop through all the images in the sweep and extract the pixels
        # from each of the images
        progress = ProgressBar(indent=4, title='Extracting pixels from sweep')
        for frame, image in enumerate(sweep):
            c, i = self._extract_image_pixels(image, frame + start)
            coords.extend(c)
            intensity.extend(i)
            progress.update(100.0 * float(frame + 1) / len(sweep))

        progress.finished('Extracted {0} strong pixels'.format(len(coords)))

        # Reuturn the pixel values
        return coords, intensity

    def _extract_image_pixels(self, image, frame):
        '''Extract pixels from a single image

        Params:
            image The image from which to extract the pixels
            frame The frame number of the current image
            trusted_range The trusted range of pixel values

        Returns:
            The list of selected pixels

        '''
        from scitbx.array_family import flex
        from dials.algorithms.peak_finding import flex_vec3_int
        import numpy
        from time import time

        # Calculate the threshold
        mask = self._threshold_strategy(image)

        # Extract the pixel indices
        index = numpy.where(mask.as_numpy_array())

        # Create the array of pixel coordinates
        z = [frame] * len(index[0])
        y = map(int, index[0])
        x = map(int, index[1])
        coords = flex_vec3_int(zip(x, y, z))

        # Get the array of pixel intensities
        intensity = flex.int(image.as_numpy_array()[index])

        # Return the pixel values
        return coords, intensity

    def _label_pixels(self, pixels, sweep):
        '''Do a connected component labelling of the pixels to get
        groups of spots.

        Params:
            pixels The pixel coordinates
            sweep The sweep object

        Returns:
            The pixel-spot mapping.

        '''
        from dials.algorithms.peak_finding import LabelPixels

        # Get the grid size needed by the labelling algorithm
        image_size = sweep.get_image_size()
        scan_size = sweep.get_array_range()[1]
        grid_size = (image_size[0], image_size[1], scan_size)

        # Label the pixels
        label_pixels = LabelPixels(grid_size)
        labels = label_pixels(pixels)

        # Return the labels
        return labels

    def _filter_spots(self, labels):
        '''Filter spots that are too small.

        Params:
            spots The input spots

        Returns:
            The filtered spots

        '''
        # Create a list of lists of pixels representing the spots
        spots = [list() for i in range(max(labels) + 1)]
        for i, l in enumerate(labels):
            spots[l].append(i)

        # Filter by spot size
        return [s for s in spots if len(s) >= self._min_spot_size]

    def _calculate_bbox(self, coords, spots, sweep):
        '''Calculate the bounding boxes for each spot.

        Params:
            coords The pixel coordinates
            spots The pixel-spot mapping
            sweep The sweep object

        Returns:
            The bounding boxes for each spot

        '''
        # Get the image dimensions
        height, width = sweep.get_image_size()
        length = sweep.get_array_range()[1]

        # Loop through all spots
        bbox = []
        for s in spots:

            # Loop through all pixels in spot and find the min/max
            # x/y/z coordinates
            b = [width, -1, height, -1, length, -1]
            for i in s:
                x, y, z = coords[i]
                if x < b[0]: b[0] = x
                if x > b[1]: b[1] = x
                if y < b[2]: b[2] = y
                if y > b[3]: b[3] = y
                if z < b[4]: b[4] = z
                if z > b[5]: b[5] = z

            # Make the bounding box and check it is valid
            b[1] += 1
            b[3] += 1
            b[5] += 1
            assert(b[0] < b[1] and b[2] < b[3] and b[4] < b[5])

            # Add to list
            bbox.append(tuple(b))

        # Return the bounding boxes
        return bbox

    def _calculate_centroids(self, coords, intensity, spots):
        '''Calculate the spot centroids.

        Params:
            coords The spot coords
            intensity The spot intensities
            spots The pixel-spot mapping

        Returns:
            (centroid position, centroid variance)

        '''
        from dials.algorithms.image.centroid import centroid_points
        from scitbx.array_family import flex

        # Initialise arrays
        centroid_pos = flex.vec3_double()
        centroid_var = flex.vec3_double()
        centroid_err = flex.vec3_double()
        centroid_tot = flex.double()

        # Loop through each spot
        for s in spots:

            # Get arrays of pixels coordinates and values
            pixel_coords = flex.vec3_double(
                [map(lambda x: x + 0.5, coords[i]) for i in s])
            pixel_values = flex.double([intensity[i] for i in s])

            # Calculate the centroid attributes
            centroid = centroid_points(pixel_values, pixel_coords)
            centroid_pos.append(centroid.mean())
            centroid_var.append(centroid.unbiased_variance())
            centroid_err.append(centroid.unbiased_standard_error_sq())
            centroid_tot.append(centroid.sum_pixels())

        # Return the centroid and variance
        return centroid_pos, centroid_var, centroid_err, centroid_tot

    def _filter_maximum_centroid(self, coords, values, spots, cpos):
        '''Filter the reflections by the distance between the maximum pixel
        value and the centroid position. If the centroid is a greater than the
        maximum separation from maximum pixel (in pixel coords) then discard.

        Params:
            coords The list of coordinates
            values The list of values
            cpos The list of centroids

        Returns:
            An index list of valid spots

        '''
        from scitbx.array_family import flex
        from scitbx import matrix
        index = []
        for si, (s, c) in enumerate(zip(spots, cpos)):
            im = flex.max_index(flex.int([values[i] for i in s]))
            xc = matrix.col(c)
            xm = matrix.col(coords[s[im]])
            if (xc - xm).length() <= self._max_separation:
                index.append(si)

        # Return the list of indices
        return index

    def _create_reflection_list(self, coords, values, spots, bbox, cpos, cvar,
                                cerr, ctot, index):
        '''Create a reflection list from the spot data.

        Params:
            coords The pixel coordinates
            values The pixel values
            spots The pixel->spot mapping
            bbox The spot bounding boxes
            cpos The centroid position
            cvar The centroid variance
            index The list of valid indices

        Returns:
            A list of reflections

        '''
        from dials.model.data import Reflection, ReflectionList
        from dials.algorithms.integration import allocate_reflection_profiles

        # Ensure the lengths are ok
        assert(len(index) > 0)
        assert(len(spots) == len(bbox))
        assert(len(spots) == len(cpos))
        assert(len(spots) == len(cvar))
        assert(len(spots) == len(cerr))
        assert(len(spots) == len(ctot))

        # Create the reflection list
        reflection_list = ReflectionList(len(index))
        for i, r in zip(index, reflection_list):
            r.bounding_box = bbox[i]
            r.centroid_position = cpos[i]
            r.centroid_variance = cerr[i]
            r.centroid_sq_width = cvar[i]
            r.intensity = ctot[i]

        # Allocate memory for the reflection profiles
        allocate_reflection_profiles(reflection_list,
            shoebox_default=0, shoebox_mask_default=0)

        # Set the pixel and mask values
        for i, r in zip(index, reflection_list):
            xs, xf, ys, yf, zs, zf = r.bounding_box
            for s in spots[i]:
                x, y, z = coords[s]
                x = x - xs
                y = y - ys
                z = z - zs
                r.shoebox[z, y, x] = values[s]
                r.shoebox_mask[z, y, x] = 1

        # Return the reflection list
        return reflection_list
