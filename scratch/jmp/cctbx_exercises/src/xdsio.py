"""A module to provide functions to read/write XDS files."""


class GxParmFile:
    """A class to read the GXPARM.XDS file used in XDS"""

    def __init__(self):
        pass

    def read_file(self, filename):
        """Read the GXPARAM.XDS file.

        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.

        :param filename: The path to the file

        """
        # Read the text from the file and split into an array of tokens
        lines = open(filename, 'r').readlines()
        tokens = [map(float, l.split()) for l in lines]

        # Read the parameters from the list of tokens
        self.starting_frame    = tokens[0][0]
        self.starting_angle    = tokens[0][1]
        self.oscillation_range = tokens[0][2]
        self.rotation_axis     = tuple(tokens[0][3:6])
        self.wavelength        = tokens[1][0]
        self.beam_vector       = tuple(tokens[1][1:4])
        self.detector_size     = tuple(map(int, tokens[2][0:2]))
        self.pixel_size        = tuple(tokens[2][2:4])
        self.detector_distance = tokens[3][0]
        self.detector_origin   = tuple(tokens[3][1:3])
        self.detector_x_axis   = tuple(tokens[4])
        self.detector_y_axis   = tuple(tokens[5])
        self.detector_normal   = tuple(tokens[6])
        self.space_group       = int(tokens[7][0])
        self.unit_cell         = tuple(tokens[7][1:7])
        self.unit_cell_a_axis  = tuple(tokens[8])
        self.unit_cell_b_axis  = tuple(tokens[9])
        self.unit_cell_c_axis  = tuple(tokens[10])


class IntegrateFile:
    """A class to read the INTEGRATE.HKL file used in XDS"""

    def __init__(self):
        """Initialise the file contents."""
        self._header = {}
        self.hkl = []
        self.iobs = []
        self.sigma = []
        self.xyzcal = []
        self.rlp = []
        self.peak = []
        self.corr = []
        self.maxc = []
        self.xyzobs = []
        self.alfbet0 = []
        self.alfbet1 = []
        self.psi = []


    def read_file(self, filename):
        """Read the INTEGRATE.HKL file.

        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.

        :param filename: The path to the file

        """
        # Read the lines from the file
        lines = open(filename, 'r').readlines()

        # Loop through the lines in the file. First off, parse the header
        # lines until we reach !END_OF_HEADER. Then parse the data lines
        # until we read !END_OF_DATA
        in_header = True
        for l in lines:
            if in_header:
                if l.strip().startswith('!END_OF_HEADER'):
                    in_header = False
                    continue
                else:
                    if not l.strip().startswith('!'):
                        continue
                    self._parse_header_line(l.strip()[1:])
            else:
                if l.strip().startswith('!END_OF_DATA'):
                    break
                else:
                    self._parse_data_line(l)

        # Set the header parameters
        self._set_header_parameters()

    def _parse_str(self, s):
        """Parse a string to either an int, float or string

        :param s: The input string
        :returns: The parsed value

        """
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return str(s)


    def _parse_value(self, value):
        """Parse the value or array of values contained in the string

        :param value: The value to parse
        :returns: The parsed value

        """
        values = value.split()
        if len(values) == 1:
            return self._parse_str(values[0])
        else:
            return tuple([self._parse_str(s) for s in values])


    def _set_header_parameters(self, ):
        """Get the parameters from the header dict

        :param name_value: The name, value parameter dict

        """
        self.space_group       = self._header['SPACE_GROUP_NUMBER']
        self.unit_cell         = self._header['UNIT_CELL_CONSTANTS']
        self.detector_size     = (self._header['NX'], self._header['NY'])
        self.pixel_size        = (self._header['QX'], self._header['QY'])
        self.starting_frame    = self._header['STARTING_FRAME']
        self.starting_angle    = self._header['STARTING_ANGLE']
        self.oscillation_range = self._header['OSCILLATION_RANGE']
        self.rotation_axis     = self._header['ROTATION_AXIS']
        self.wavelength        = self._header['X-RAY_WAVELENGTH']
        self.beam_vector       = self._header['INCIDENT_BEAM_DIRECTION']
        self.detector_x_axis   = self._header['DIRECTION_OF_DETECTOR_X-AXIS']
        self.detector_x_axis   = self._header['DIRECTION_OF_DETECTOR_Y-AXIS']
        self.detector_origin   = (self._header['ORGX'], self._header['ORGY'])
        self.detector_distance = self._header['DETECTOR_DISTANCE']
        self.unit_cell_a_axis  = self._header['UNIT_CELL_A-AXIS']
        self.unit_cell_b_axis  = self._header['UNIT_CELL_B-AXIS']
        self.unit_cell_c_axis  = self._header['UNIT_CELL_C-AXIS']
        self._header = None


    def _parse_header_line(self, line):
        """Parse a line that has been identified as a header line

        :param line: The line to parse

        """
        name_value = line.split('=')
        if (len(name_value) < 2):
            return

        name = name_value[0]
        if (len(name_value) > 2):
            for i in range(1, len(name_value)-1):
                value_name = name_value[i].split()
                value = ''.join(value_name[:-1])
                self._header[name] = self._parse_value(value)
                name = value_name[-1]

        value = name_value[-1]
        self._header[name] = self._parse_value(value)


    def _parse_data_line(self, line):
        """Parse a data line from the Integrate.hkl file

        :param line: The line to parse

        """
        # Split the tokens
        tokens = line.split()
        tokens = map(int, tokens[0:3]) + map(float, tokens[3:])

        # Get the reflection information and append to the lists
        self.hkl    .append(tuple(tokens[0:3]))
        self.iobs   .append(tokens[3])
        self.sigma  .append(tokens[4])
        self.xyzcal .append(tuple(tokens[5:8]))
        self.rlp    .append(tokens[8])
        self.peak   .append(tokens[9])
        self.corr   .append(tokens[10])
        self.maxc   .append(tokens[11])
        self.xyzobs .append(tuple(tokens[12:15]))
        self.alfbet0.append(tuple(tokens[15:17]))
        self.alfbet1.append(tuple(tokens[17:19]))
        self.psi    .append(tokens[19])
