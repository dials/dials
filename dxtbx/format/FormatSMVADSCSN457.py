#!/usr/bin/env python
# FormatSMVADSCSN457.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for example on Australian Synchrotron SN 457
# which has reversed phi.

import time

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVADSCSN457(FormatSMVADSC):
    '''A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 457.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this is ADSC SN 457.'''

        if FormatSMVADSC.understand(image_file) == 0:
            return 0

        # check this is detector serial number 457

        size, header = FormatSMVADSC.get_smv_header(image_file)

        if int(header['DETECTOR_SN']) != 457:
            return 0

        return 3

    def __init__(self, image_file):
        '''Initialise the image structure from the given file, including a
        proper model of the experiment.'''

        assert(FormatSMVADSC.understand(image_file) > 0)

        FormatSMV.__init__(self, image_file)

        return

    def _goniometer(self):
        '''Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header.'''

        return self._goniometer_factory.known_axis((-1, 0, 0))

    # FIXME surely I don't need the code which follows which just reproduces
    # standard ADSC model?

    # FIXME find a test case, remove this.

    def _detector(self):
        '''Return a model for a simple detector, presuming no one has
        one of these on a two-theta stage. Assert that the beam centre is
        provided in the Mosflm coordinate frame.'''

        distance = float(self._header_dictionary['DISTANCE'])
        beam_x = float(self._header_dictionary['BEAM_CENTER_X'])
        beam_y = float(self._header_dictionary['BEAM_CENTER_Y'])
        pixel_size = float(self._header_dictionary['PIXEL_SIZE'])
        image_size = (float(self._header_dictionary['SIZE1']),
                      float(self._header_dictionary['SIZE2']))
        overload = 65535
        underload = 0

        return self._detector_factory.simple(
            'CCD', distance, (beam_y, beam_x), '+x', '-y',
            (pixel_size, pixel_size), image_size, (underload, overload), [])

    def _beam(self):
        '''Return a simple model for the beam.'''

        wavelength = float(self._header_dictionary['WAVELENGTH'])

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        '''Return the scan information for this image.'''

        format = self._scan_factory.format('SMV')
        exposure_time = float(self._header_dictionary['TIME'])
        epoch =  time.mktime(time.strptime(self._header_dictionary['DATE']))
        osc_start = float(self._header_dictionary['OSC_START'])
        osc_range = float(self._header_dictionary['OSC_RANGE'])

        return self._scan_factory.single(
            self._image_file, format, exposure_time,
            osc_start, osc_range, epoch)

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatSMVADSC.understand(arg)
