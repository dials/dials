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
        
    def get_goniometer(self):
        """Get the goniometer parameters from the file
        
        Returns:
            An instance of the goniometer struct
        
        """
        from dials.equipment import Goniometer
        from scitbx import matrix
        #from math import pi
        #d2r = pi / 180.0
        return Goniometer(matrix.col(self.rotation_axis).normalize().elems,
                          self.starting_angle,# * d2r,
                          self.oscillation_range,# * d2r,
                          int(self.starting_frame))
        
    def get_beam(self):
        """Get the beam parameters from the file
        
        Returns:
            An instance of the beam struct
        
        """
        from dials.equipment import Beam
        from scitbx import matrix
        return Beam(matrix.col(self.beam_vector).normalize() / self.wavelength, 
                    self.wavelength)
        
    def get_detector(self):
        """Get the detector parameters from the file
        
        Returns:
            An instance of the detector struct
        
        """    
        from dials.equipment import Detector
        from scitbx import matrix
        return Detector(matrix.col(self.detector_x_axis).normalize().elems,
                        matrix.col(self.detector_y_axis).normalize().elems,
                        matrix.col(self.detector_normal).normalize().elems,
                        self.detector_origin,
                        self.pixel_size,
                        self.detector_size,
                        self.detector_distance)
                                          
    def get_ub_matrix(self):
        """Get the UB matrix
        
        Returns:
            The UB matrix
        
        """
        from scitbx import matrix
        return matrix.sqr(self.unit_cell_a_axis + 
                          self.unit_cell_b_axis + 
                          self.unit_cell_c_axis)
                          
    def get_unit_cell(self):
        from cctbx import uctbx
        return uctbx.unit_cell(orthogonalization_matrix = self.get_ub_matrix())

    def get_space_group_type(self):
        from cctbx import sgtbx
        return sgtbx.space_group_type(sgtbx.space_group(
                            sgtbx.space_group_symbols(self.space_group).hall()))
        
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


    def _set_header_parameters(self):
        """Get the parameters from the header dict
        
        :param name_value: The name, value parameter dict
        
        """
        from scitbx import matrix
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
        self.detector_y_axis   = self._header['DIRECTION_OF_DETECTOR_Y-AXIS']
        self.detector_origin   = (self._header['ORGX'], self._header['ORGY'])
        self.detector_distance = self._header['DETECTOR_DISTANCE']
        self.unit_cell_a_axis  = self._header['UNIT_CELL_A-AXIS']
        self.unit_cell_b_axis  = self._header['UNIT_CELL_B-AXIS']
        self.unit_cell_c_axis  = self._header['UNIT_CELL_C-AXIS']
        self.sigma_divergence  = self._header['BEAM_DIVERGENCE_E.S.D.']
        self.sigma_mosaicity   = self._header['REFLECTING_RANGE_E.S.D.']
 
        # Normalize a few vectors
        self.detector_x_axis   = tuple(matrix.col(self.detector_x_axis).normalize())
        self.detector_y_axis   = tuple(matrix.col(self.detector_y_axis).normalize())
        self.detector_normal   = tuple(matrix.col(self.detector_x_axis).cross(
                                    matrix.col(self.detector_y_axis)))
        self.rotation_axis     = tuple(matrix.col(self.rotation_axis).normalize())
        self.beam_vector       = tuple(matrix.col(self.beam_vector).normalize() / self.wavelength)
        del(self._header)


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
        
    def get_goniometer(self):
        """Get the goniometer parameters from the file
        
        Returns:
            An instance of the goniometer struct
        
        """
        from dials.equipment import Goniometer
        #from math import pi
        #d2r = pi / 180.0
        return Goniometer(self.rotation_axis,
                          self.starting_angle,
                          self.oscillation_range,
                          int(self.starting_frame))
        
    def get_beam(self):
        """Get the beam parameters from the file
        
        Returns:
            An instance of the beam struct
        
        """
        from dials.equipment import Beam
        from scitbx import matrix
        return Beam(matrix.col(self.beam_vector).normalize() / self.wavelength, 
                    self.wavelength)
        
    def get_detector(self):
        """Get the detector parameters from the file
        
        Returns:
            An instance of the detector struct
        
        """    
        from dials.equipment import Detector
        return Detector(self.detector_x_axis,
                        self.detector_y_axis,
                        self.detector_normal,
                        self.detector_origin,
                        self.pixel_size,
                        self.detector_size,
                        self.detector_distance)
                                          
    def get_ub_matrix(self):
        """Get the UB matrix
        
        Returns:
            The UB matrix
        
        """
        from scitbx import matrix
        return matrix.sqr(self.unit_cell_a_axis + 
                          self.unit_cell_b_axis + 
                          self.unit_cell_c_axis)
  
class XYCorrection:
    
    def __init__(self):
        pass
        
    def read_file(self, filename):
        """Read the CBF correction file"""
        import pycbf
        self.cbf_handle = pycbf.cbf_handle_struct()
        self.cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
        self.cbf_handle.rewind_datablock()

    def get_correction_array(self):
        """Get the correction array from the file"""
        import numpy
        
        # Select the first datablock and rewind all the categories
        self.cbf_handle.select_datablock(0)
        self.cbf_handle.select_category(0)
        self.cbf_handle.select_column(2)
        self.cbf_handle.select_row(0)
    
        # Check the type of the element to ensure it's a binary
        # otherwise raise an exception
        type = self.cbf_handle.get_typeofvalue()
        if type.find('bnry') > -1:

            # Read the image data into an array
            image_string = self.cbf_handle.get_integerarray_as_string()
            image = numpy.fromstring(image_string, numpy.int32)
            
            # Get the array parameters
            parameters = self.cbf_handle.get_integerarrayparameters_wdims()
            image_size = (parameters[10], parameters[9])

            # Resize the image
            image.shape = (image_size)

        else:
            raise TypeError('Can\'t find image')        

        # Return the image
        return image

    def get_correction(self, dim):
        """Get the correction at each pixel."""
        import numpy
        
        # Get the raw array
        raw_array = self.get_correction_array()

        # Ensure out dimensions are ok
        if raw_array.shape[0] * 4 < dim[0] or raw_array.shape[1] * 4 < dim[1]:
            raise ValueError("Dimensions are incompatible")

        # Create the array of the given dimension
        correction = numpy.zeros(dim, dtype=numpy.float64)

        # Loop through all pixels and get the correction
        i1 = numpy.array([range(dim[1])] * dim[0], dtype=numpy.int32)
        j1 = numpy.array([range(dim[0])] * dim[1], dtype=numpy.int32).transpose()
        i2 = numpy.divide(i1, 4)
        j2 = numpy.divide(j1, 4)
        correction[j1,i1] = raw_array[j2,i2] / 10.0
                
        # Return the array of corrections
        return correction
        
