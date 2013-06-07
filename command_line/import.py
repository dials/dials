#!/usr/bin/env python
#
# import.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class RawSweepImporter(object):
    ''' Class to import sweep from raw filenames. '''

    def __init__(self, filenames):
        ''' Set the filenames. '''
        self.filenames = filenames

    def __call__(self):
        ''' Get the imageset data.'''
        from dxtbx.imageset import ImageSetFactory, ImageSet, ImageSweep

        # Get the imagsets from the input files
        imagesets = ImageSetFactory.new(self.filenames)

        # Ensure we have 1 imageset
        if len(imagesets) != 1:
            raise RuntimeError(
                'Requires exactly 1 sweep, got {0}'.format(len(imagesets)))

        # Ensure its a sweep
        if not isinstance(imagesets[0], ImageSweep):
            raise RuntimeError('Dials requires an ImageSweep nor an ImageSet')

        # Return the sweep
        return imagesets[0]


class XdsSweepImporter(object):

    def __init__(self, directory):
        pass

    def __call__(self):
        return None


class XdsCrystImporter(object):
    ''' Class to import crystal from XDS files.'''

    def __init__(self, directory, file_choice):
        '''Init the class.'''

        self.directory = directory
        self.file_choice = file_choice

    def get_filename(self, choice):
        '''Get filename from choice.'''
        import os

        # Return None for bad choice
        if choice == None:
            return None

        # Get path for choice
        filename = {
          'XPARM' : 'XPARM.XDS',
          'GXPARM' : 'GXPARM.XDS',
          'INTEGRATE' : 'INTEGRATE.HKL',
          'XDS_ASCII' : 'XDS_ASCII.HKL'
        }
        return os.path.join(self.directory, filename[choice])

    def find_best_file(self):
        ''' Find the best available file.'''
        import os
        if os.path.exists(self.get_filename('XDS_ASCII')):
            return self.get_filename('XDS_ASCII')
        elif os.path.exists(self.get_filename('INTEGRATE.HKL')):
            return self.get_filename('INTEGRATE.HKL')
        elif os.path.exists(self.get_filename('GXPARM.XDS')):
            return self.get_filename('GXPARM.XDS')
        elif os.path.exists(self.get_filename('XPARM')):
            return self.get_filename('XPARM')
        else:
            raise RuntimeError('No XDS files found.')

    def __call__(self):
        ''' Import the crystal. '''
        import os
        from dials.model.serialize import xds

        # Get the filename if given a choice
        filename = self.get_filename(self.file_choice)

        # If the choice doesn't exist then find best file
        if filename == None or not os.path.exists(filename):
            filename = self.find_best_file()

        # Return the crystal model
        return xds.to_crystal(filename)


class Script(ScriptRunner):
    '''Class to run the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # Specify the command line options
        usage = "usage: %prog [options] [/path/to/image/files]"

        ScriptRunner.__init__(self,
            usage=usage,
            show_config_option=False,
            save_config_option=False)

        # Add options to do with XDS
        self.config().add_option(
            '--xds',
            dest = 'xds_dir',
            type = 'string', default = None,
            help = 'Directory containing XDS files.')

        # XDS crystal choice
        self.config().add_option(
            '--xds-crystal',
            dest = 'xds_crystal',
            type = 'choice',
            choices = ['XPARM', 'GXPARM', 'INTEGRATE', 'XDS_ASCII'],
            default = None,
            help = 'Choose the XDS file to use for the Crystal model')

    def main(self, params, options, args):
        ''' Run the script.'''

        from dials.model.serialize import dump
        from dials.util.options import ConfigWriter

        # Add raw importers
        if len(args) > 0:
            import_sweep = RawSweepImporter(args)
            import_cryst = lambda: None
            import_param = lambda: params
        elif options.xds_dir:
            import_sweep = XdsSweepImporter(options.xds_dir)
            import_cryst = XdsCrystImporter(
                options.xds_dir, options.xds_crystal)
            import_param = lambda: params
        else:
            self.config().print_help()
            return

        # Set the output filenames
        sweep_filename = 'sweep.json'
        cryst_filename = 'crystal.json'
        param_filename = 'params.phil'

        # Import the sweep, crystal, params from the input data representation
        sweep = import_sweep()
        cryst = import_cryst()
        param = import_param()

        # Serialize the sweep, crystal and params
        if sweep:
            dump.sweep(sweep, sweep_filename)
            print 'Saved sweep to {0}'.format(sweep_filename)
        if cryst:
            dump.crystal(cryst, cryst_filename)
            print 'Saved crystal to {0}'.format(cryst_filename)
        if param:
            writer = ConfigWriter(self.config().system_phil())
            writer.write(param, param_filename)
            print 'Saved parameters to {0}'.format(param_filename)


if __name__ == '__main__':
    script = Script()
    script.run()
