#!/usr/bin/env python
#
# extract_spot_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.


class ScriptRunner(object):
    '''Class to run script.'''

    def __init__(self, spot_filename, xparm_filename,
                 output_filename, include_invalid):
        '''Setup the script.'''

        # Filename data
        self.spot_filename = spot_filename
        self.xparm_filename = xparm_filename
        self.output_filename = output_filename
        self.include_invalid = include_invalid

    def __call__(self):
        '''Run the script.'''
        from iotbx.xds import spot_xds
        from dials.model.data import Reflection, ReflectionList
        from dials.util.command_line import Command
        import dxtbx

        # Read the SPOT.XDS file
        Command.start('Reading SPOT.XDS')
        handle = spot_xds.reader()
        handle.read_file(self.spot_filename)

        # Get the info
        centroid = handle.centroid
        intensity = handle.intensity
        try:
            miller_index = handle.miller_index
        except AttributeError:
            miller_index = None

        Command.end('Read {0} spots from SPOT.XDS file.'.format(len(centroid)))

        # Read the models from the xparm file
        Command.start('Reading models from {0}'.format(self.xparm_filename))
        models = dxtbx.load(self.xparm_filename)
        self.detector = models.get_detector()
        self.scan = models.get_scan()
        self.pixel_size = self.detector.get_pixel_size()
        self.oscillation_range = self.scan.get_oscillation(deg=False)
        Command.end('Read models from {0}'.format(self.xparm_filename))

        # Create the reflection list
        Command.start('Creating reflection list')
        if miller_index:
            rlist = ReflectionList()
            for c, i, h in zip(centroid, intensity, miller_index):
                if self.include_invalid == True or tuple(h) != (0, 0, 0):
                    r = self._get_reflection_data(c, i, h)
                    rlist.append(r)

        else:
            rlist = ReflectionList()
            for r, c, i in zip(rlist, centroid, intensity):
                r = self._get_reflection_data(c, i, (0, 0, 0))
                rlist.append(r)

        Command.end('Created reflection list')

        # Save the reflection list
        if self.output_filename != None:
            import pickle
            Command.start('Saving reflections to {0}'.format(
                self.output_filename))
            pickle.dump(rlist, open(self.output_filename, 'wb'))
            Command.end('Saved reflections to {0}'.format(
                self.output_filename))

    def _get_reflection_data(self, pos_px, intensity, hkl):
        '''Create the reflection data'''
        from dials.model.data import Reflection
        from dials.algorithms.centroid import centroid_px_to_mm

        r = Reflection()

        # Centroid position/variance in pixels
        sqw_px = (1.0, 1.0, 1.0)
        var_px = (1.0 + sqw_px[0] / intensity,
                  1.0 + sqw_px[1] / intensity,
                  1.0 + sqw_px[2] / intensity)

        # Get the centroid in mm/rad
        pos_mm, var_mm, sqw_mm = centroid_px_to_mm(self.detector,
            self.scan, pos_px, var_px, sqw_px)

        # Put everything into the reflection struct
        r.centroid_position = pos_mm
        r.centroid_variance = var_mm
        r.centroid_sq_width = sqw_mm
        r.rotation_angle    = pos_mm[2]
        r.image_coord_mm    = (pos_mm[0], pos_mm[1])
        r.image_coord_px    = (pos_px[0], pos_px[1])
        r.frame_number      = pos_px[2]
        r.intensity         = intensity
        r.miller_index      = hkl

        # Return reflection
        return r


if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage  = "usage: %prog [options] "
    usage += "/path/to/SPOT.XDS "
    usage += "/path/to/XPARM.XDS"

    # Create an option parser
    parser = OptionParser(usage)

    # Add command line options
    parser.add_option('-z', '--zero-hkl',
                      dest='zero_hkl', action="store_true", default=False,
                      help='Include invalid reflections')
    parser.add_option('-o', '--output-file',
                      dest='output_file', type="string", default=None,
                      help='Destination filename for reflections')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 2:
        parser.print_help()
    else:
        # Initialise the script runner
        runner = ScriptRunner(
            spot_filename=args[0],
            xparm_filename=args[1],
            output_filename=options.output_file,
            include_invalid=options.zero_hkl)

        # Run the script
        runner()
