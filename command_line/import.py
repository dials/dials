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


class XDSFile(object):
    '''Class to find XDS config files.'''

    @staticmethod
    def get_filename(directory, choice):
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
        return os.path.join(directory, filename[choice])

    @staticmethod
    def find_best_file(directory):
        ''' Find the best available file.'''
        import os
        if os.path.exists(XDSFile.get_filename(directory, 'XDS_ASCII')):
            return XDSFile.get_filename(directory, 'XDS_ASCII')
        elif os.path.exists(XDSFile.get_filename(directory, 'INTEGRATE')):
            return XDSFile.get_filename(directory, 'INTEGRATE')
        elif os.path.exists(XDSFile.get_filename(directory, 'GXPARM')):
            return XDSFile.get_filename(directory, 'GXPARM')
        elif os.path.exists(XDSFile.get_filename(directory, 'XPARM')):
            return XDSFile.get_filename(directory, 'XPARM')
        else:
            return None


class XdsSweepImporter(object):
    ''' Class to import sweep from XDS files.'''

    def __init__(self, directory, file_choice):
        '''Init the class.'''

        self.directory = directory
        self.file_choice = file_choice

    def __call__(self):
        ''' Import the sweep. '''
        import os
        from dials.model.serialize import xds

        # Set the input filename
        input_filename = os.path.join(self.directory, 'XDS.INP')

        # Get the filename if given a choice
        extra_filename = XDSFile.get_filename(self.directory, self.file_choice)

        # If the choice doesn't exist then find best file
        if extra_filename == None or not os.path.exists(extra_filename):
            extra_filename = XDSFile.find_best_file(self.directory)

        # Return the crystal model
        return xds.to_sweep(input_filename, extra_filename)


class XdsCrystImporter(object):
    ''' Class to import crystal from XDS files.'''

    def __init__(self, directory, file_choice):
        '''Init the class.'''

        self.directory = directory
        self.file_choice = file_choice

    def __call__(self):
        ''' Import the crystal. '''
        import os
        from dials.model.serialize import xds

        # Get the filename if given a choice
        filename = XDSFile.get_filename(self.directory, self.file_choice)

        # If the choice doesn't exist then find best file
        if filename == None or not os.path.exists(filename):
            filename = XDSFile.find_best_file(self.directory)
            if filename == None:
                raise RuntimeError('Unable to find crystal file.')

        # Return the crystal model
        return xds.to_crystal(filename)


class XdsParamImporter(object):
    ''' Class to import parameters from XDS files. '''

    def __init__(self, directory):
        ''' Init the class. '''
        self.directory = directory

    def __call__(self):
        ''' Import the parameters '''
        return ''


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

        # XDS sweep choice
        self.config().add_option(
            '--xds-sweep',
            dest = 'xds_sweep',
            type = 'choice',
            choices = ['XPARM', 'GXPARM', 'INTEGRATE', 'XDS_ASCII'],
            default = None,
            help = 'Choose the XDS file to use for the geometry models: '
                   'XPARM|GXPARM|INTEGRATE|XDS_ASCII')

        # XDS crystal choice
        self.config().add_option(
            '--xds-crystal',
            dest = 'xds_crystal',
            type = 'choice',
            choices = ['XPARM', 'GXPARM', 'INTEGRATE', 'XDS_ASCII'],
            default = None,
            help = 'Choose the XDS file to use for the Crystal model: '
                   'XPARM|GXPARM|INTEGRATE|XDS_ASCII')

        # Add options to do with XDS
        self.config().add_option(
            '--sweep-filename',
            dest = 'sweep_filename',
            type = 'string', default = "sweep.json",
            help = 'Filename for output sweep file.')

        # Add options to do with XDS
        self.config().add_option(
            '--crystal-filename',
            dest = 'crystal_filename',
            type = 'string', default = "crystal.json",
            help = 'Filename for output crystal file.')

        # Add options to do with XDS
        self.config().add_option(
            '--param-filename',
            dest = 'param_filename',
            type = 'string', default = "param.phil",
            help = 'Filename for output param file.')

    def main(self, params, options, args):
        ''' Run the script.'''

        from dials.model.serialize import dump
        from dials.util.options import ConfigWriter

        # No source!
        sweep_source = "nowhere!"
        cryst_source = "nowhere!"
        param_source = "system parameters"

        # Add raw importers
        if len(args) > 0:

            # Get sweep from image file input
            import_sweep = RawSweepImporter(args)
            import_cryst = lambda: None
            import_param = lambda: params
            sweep_source = "image file data"

            # See if we can get crystal from XDS
            if options.xds_dir:
                import_cryst = XdsCrystImporter(
                    options.xds_dir,
                    options.xds_crystal)
                import_param = lambda: params
                cryst_source = "xds configuration files"

        # Get sweep and crystal from XDS
        elif options.xds_dir:
            import_sweep = XdsSweepImporter(
                options.xds_dir,
                options.xds_sweep)
            import_cryst = XdsCrystImporter(
                options.xds_dir,
                options.xds_crystal)
            import_param = XdsParamImporter(
                options.xds_dir)
            sweep_source = "xds configuration files"
            cryst_source = "xds configuration files"
            param_source = "xds configuration files"
        else:
            self.config().print_help()
            return

        # Print some information
        print "Importing data from the following sources:"
        print " - Sweep from {0}".format(sweep_source)
        print " - Crystal from {0}".format(cryst_source)
        print " - Parameters from {0}".format(param_source)

        # Import the sweep, crystal, params from the input data representation
        sweep = import_sweep()
        cryst = import_cryst()
        param = import_param()

        # Serialize the sweep, crystal and params
        if sweep:
            dump.sweep(sweep, options.sweep_filename)
            print 'Saved sweep to {0}'.format(options.sweep_filename)
        if cryst:
            dump.crystal(cryst, options.crystal_filename)
            print 'Saved crystal to {0}'.format(options.crystal_filename)
        if param:
            writer = ConfigWriter(self.config().system_phil())
            writer.write(param, options.param_filename)
            print 'Saved parameters to {0}'.format(options.param_filename)


if __name__ == '__main__':
    script = Script()
    script.run()
