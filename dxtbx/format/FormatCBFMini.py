#!/usr/bin/env python
# FormatCBFMini.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of miniCBF format - as used with Dectris detectors -
# this will read the header and populate a dictionary of the keyword / value
# pairs.

from dxtbx.format.FormatCBF import FormatCBF

class FormatCBFMini(FormatCBF):
    '''An image reading class for mini CBF format images i.e. those from
    Dectris, which will read the header into a dictionary.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an CBF format image, i.e. we can
        make sense of it. N.B. in situations where there is both a full and
        mini CBF header this will return 0 i.e. deprecating itself in favour
        of the full CBF reader.'''

        if FormatCBF.understand(image_file) == 0:
            return 0

        header = FormatCBF.get_cbf_header(image_file)

        if '_diffrn.id' in header and '_diffrn_source' in header:
            return 0

        for record in header.split('\n'):
            if '_array_data.header_convention' in record and \
                   'PILATUS' in record:
                return 2
            if '_array_data.header_convention' in record and \
                   'SLS' in record:
                return 2
            if '_array_data.header_convention' in record and \
                   '?' in record:
                return 2

        return 0

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(FormatCBFMini.understand(image_file) > 0)

        FormatCBF.__init__(self, image_file)

        return

    def _start(self):
        '''Open the image file, read the image header, copy it into a
        dictionary for future reference.'''

        FormatCBF._start(self)

        cif_header = FormatCBF.get_cbf_header(self._image_file)

        self._cif_header_dictionary = { }

        for record in cif_header.split('\n'):
            if not '#' in record[:1]:
                continue

            if len(record[1:].split()) <= 2 and record.count(':') == 2:
                self._cif_header_dictionary['timestamp'] = record[1:].strip()
                continue

            tokens = record.replace('=', '').replace(':', '').split()[1:]

            self._cif_header_dictionary[tokens[0]] = ' '.join(tokens[1:])

        for record in self._mime_header.split('\n'):
            if not record.strip():
                continue
            token, value = record.split(':')
            self._cif_header_dictionary[token.strip()] = value.strip()

        return
