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


class ExtractSpots(object):
    ''' Class to find spots in an image and extract them into shoeboxes. '''

    def __init__(self, threshold_image):
        ''' Initialise the class with the strategy

        Params:
            threshold_image The image thresholding strategy

        '''
        # Set the required strategies
        self.threshold_image = threshold_image

    def __call__(self, sweep):
        ''' Find the spots in the sweep

        Params:
            sweep The sweep to process

        Returns:
            The list of spot shoeboxes

        '''
        from dials.util.command_line import ProgressBar, Command
        from dials.algorithms.image.connected_components import LabelImageStack
        from dials.array_family import flex

        # Construct the pixel labeller
        label = [LabelImageStack(panel.get_image_size()[::-1])
            for panel in sweep.get_detector()]

        # Loop through all the images in the sweep and extract the pixels
        # from each of the images
        progress = ProgressBar(title='Extracting pixels from sweep')
        for frame, image in enumerate(sweep):

            # Ensure image is a tuple of images
            if not isinstance(image, tuple):
                image = (image,)

            # Create the mask by thresholding the image. Then add the mask
            # and the image to the pixel labeller
            for lb, im in zip(label, image):
                lb.add_image(im, self.threshold_image(im))

            # Update progress bar
            progress.update(100.0 * float(frame + 1) / len(sweep))

        # Finish the progess bar
        progress.finished('Extracted {0} strong pixels'.format(
          sum([len(lb.values()) for lb in label])))

        # Extract the shoeboxes
        Command.start('Extracting spots from pixels')
        shoeboxes = flex.shoebox()
        for i, lb in enumerate(label):
            shoeboxes.extend(flex.shoebox(lb, i))
        Command.end('Extracted {0} spots from pixels'.format(len(shoeboxes)))

        # Return the shoeboxes
        return shoeboxes


class SpotFinder(object):
    ''' A class to do spot finding and filtering. '''

    def __init__(self, find_spots=None, filter_spots=None, scan_range=None):
        ''' Initialise the class. '''

        # Set the spot finding and filter functions
        assert(find_spots != None and filter_spots != None)
        self.find_spots = find_spots
        self.filter_spots = filter_spots

        # Set the scan range
        self.scan_range = scan_range

    def __call__(self, sweep):
        ''' Do the spot finding '''
        from dials.model.data import ReflectionList
        from dials.array_family import flex
        from dials.util.command_line import Command

        # Get list of scan ranges
        if not self.scan_range:
            scan_range = [(0, len(sweep))]
        else:
            scan_range = self.scan_range

        # Get spots from bits of scan
        spots_all = []
        for scan in scan_range:
            j0, j1 = scan
            assert(j1 > j0 and j0 >= 0 and j1 <= len(sweep))
            print '\nFinding spots in image {0} to {1}...'.format(j0, j1)

            # Find the spots
            spots = self.find_spots(sweep[j0:j1])

            # Add the spots to the list
            spots_all.extend(spots)

        # Get the list of shoeboxes
        shoeboxes = flex.shoebox(spots_all)

        # Calculate the spot centroids
        Command.start('Calculating {0} spot centroids'.format(len(shoeboxes)))
        centroid = shoeboxes.centroid_valid();
        Command.end('Calculated {0} spot centroids'.format(len(shoeboxes)))

        # Calculate the spot intensities
        Command.start('Calculating {0} spot intensities'.format(len(shoeboxes)))
        intensity = shoeboxes.summed_intensity_valid();
        Command.end('Calculated {0} spot intensities'.format(len(shoeboxes)))

        # Create the observations
        observed = flex.observation(shoeboxes.panels(), centroid, intensity)

        # Filter the reflections and select only the desired spots
        flags = self.filter_spots(None,
            sweep=sweep,
            observations=observed,
            shoeboxes=shoeboxes)
        observed = observed.select(flags)
        shoeboxes = shoeboxes.select(flags)

        # Return as a reflection list
        return ReflectionList(observed, shoeboxes)
