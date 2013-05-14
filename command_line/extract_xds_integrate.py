#!/usr/bin/env python
#
# extract_xds_integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, xparm_filename, integrate_filename,
                 image_filenames, output_filename):
        '''Setup the script.'''

        # Filename data
        self.xparm_filename = xparm_filename
        self.integrate_filename = integrate_filename
        self.image_filenames = image_filenames
        self.output_filename = output_filename

    def __call__(self):
        '''Run the script.'''
        from dials.algorithms.integration import BBoxCalculator
        from dials.algorithms.integration import extract_reflection_profiles
        from iotbx.xds import integrate_hkl
        from dials.model.data import Reflection, ReflectionList
        from dials.util.command_line import Command
        from dxtbx.imageset import ImageSetFactory
        import dxtbx
        from math import pi
        from scitbx import matrix

        # Read the SPOT.XDS file
        Command.start('Reading INTEGRATE.HKL')
        handle = integrate_hkl.reader()
        handle.read_file(self.integrate_filename)
        hkl    = handle.hkl
        xyzcal = handle.xyzcal
        xyzobs = handle.xyzobs
        sigma_divergence = handle.sigma_divergence
        sigma_mosaicity = handle.sigma_mosaicity
        Command.end('Read {0} reflections from INTEGRATE.HKL file.'.format(
            len(hkl)))

        # Set the number of frames
        num_frames = len(self.image_filenames)

        # Read the models from the input file
        print "Reading: \"{0}\"".format(self.xparm_filename)
        models = dxtbx.load(self.xparm_filename)
        beam = models.get_beam()
        detector = models.get_detector()
        gonio = models.get_goniometer()
        scan = models.get_scan()
        first_image = scan.get_image_range()[0]
        image_range = (first_image, first_image + num_frames)
        scan.set_image_range(image_range)

        # Calculate the bbox sizes
        n_sigma = 5.0
        delta_divergence = n_sigma * sigma_divergence * pi / 180.0
        delta_mosaicity = n_sigma * sigma_mosaicity * pi / 180.0

        # Create the reflection list
        rlist = ReflectionList(len(hkl))
        for r, h, x1, x2 in zip(rlist, hkl, xyzcal, xyzobs):
            r.miller_index = h
            r.image_coord_px = x1[0:2]
            r.frame_number = x1[2]
            r.centroid_position = x2

            # Calculate the beam vectors
            s1 = matrix.col(detector.get_pixel_lab_coord(r.image_coord_px))
            s1 /= beam.get_wavelength()
            r.beam_vector = s1

            # Calculate rotation anlge
            r.rotation_angle = scan.get_angle_from_array_index(
                r.frame_number, deg=False)

        # Calculate the bounding box for each reflection
        calculate_bbox = BBoxCalculator(beam, detector, gonio, scan,
            delta_divergence, delta_mosaicity)

        # Calculate the bounding bxoes of all the reflections
        Command.start('Calculating bounding boxes')
        calculate_bbox(rlist)
        Command.end('Calculated {0} bounding boxes.'.format(len(rlist)))

        # Extract the reflection profiles
        Command.start('Extracting reflection profiles')
        sweep = ImageSetFactory.new(self.image_filenames)[0]
        reflections = extract_reflection_profiles(sweep, rlist)
        Command.end('Extracted {0} reflection profiles.'.format(len(rlist)))

        # Filter valid reflections
        Command.start('Filtering reflections')
        rlist = self.filter_reflection_profiles(rlist, detector, scan)
        Command.end('Filtered {0} reflections'.format(len(rlist)))

        # Save the reflection list
        if self.output_filename != None:
            import pickle
            Command.start('Saving reflections to {0}'.format(
                self.output_filename))
            pickle.dump(rlist, open(self.output_filename, 'wb'))
            Command.end('Saved reflections to {0}'.format(
                self.output_filename))

    def filter_reflection_profiles(self, reflections, detector, scan):
        '''Filter invalid reflections.'''
        from dials.model.data import ReflectionList
        width, height = detector.get_image_size()
        nframes = scan.get_array_range()[1]
        valid = []
        for r in reflections:
            bbox = r.bounding_box
            if bbox[0] < 0 or bbox[2] < 0 or bbox[4] < 0:
                continue
            if bbox[1] > width or bbox[3] > height or bbox[5] > nframes:
                continue
            shoebox = r.shoebox
            if not (shoebox >= 0).all_eq(True):
                continue
            valid.append(r)
        return ReflectionList(valid)

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] "
    usage += "/path/to/XPARM.XDS "
    usage += "/path/to/INTEGRATE.HKL "
    usage += "/path/to/image/files.cbf"

    # Create an option parser
    parser = OptionParser(usage)

    # Add command line options
    parser.add_option('-o', '--output-file',
                      dest='output_file', type="string", default=None,
                      help='Destination filename for reflections')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 3:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(
            xparm_filename=args[0],
            integrate_filename=args[1],
            image_filenames=args[2:],
            output_filename=options.output_file)

        # Run the script
        runner()
