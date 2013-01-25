
class SpotPredictor(object):
    """A class to perform spot prediction
    
    Usage:
    
        spot_predictor = SpotPredictor(beam, goniometer, ub_matrix, d_min,
                                       unit_cell, space_group_type)
                                       
        spot_predictor.predict()
        hkl = spot_predictor.miller_indices
        phi = spot_predictor.rotation_angles
        s1  = spot_predictor.beam_vectors
        xyz = spot_predictor.image_volume_coordinates
    
    """

    def __init__(self, beam, goniometer, detector, ub_matrix, d_min, unit_cell, 
                 space_group_type):
        """Initialise the spot predictor.
        
        Args:
            beam The beam structure
            goniometer The goniometer structure
            detector The detector structure
            ub_matrix The UB matrix
            d_min The resolution
            unit_cell The unit_cell structure
            space_group_type The space group structure
        
        """
        # Reset the internal arrays
        self.reset()

        # Set the equipment
        self._beam = beam
        self._goniometer = goniometer
        self._detector = detector

        # Set the crystal stuff
        self._unit_cell = unit_cell
        self._space_group_type = space_group_type
        self._d_min = d_min
        self._ub_matrix = ub_matrix.inverse()

    @property
    def miller_indices(self):
        """Get the calculated miller indices."""
        return self._miller_indices
    
    @property
    def rotation_angles(self):
        """Get the calculated rotation angles."""
        return self._rotation_angles
        
    @property
    def beam_vectors(self):
        """Get the calculated beam vectors."""
        return self._beam_vectors
        
    @property
    def image_volume_coordinates(self):
        """Get the calculated image volume coordinates."""
        return self._image_coords
    
    def reset(self):
        """Reset the internal arrays to None."""
        self._miller_indices = None
        self._rotation_angles = None
        self._beam_vectors = None
        self._image_coords = None
    
    def predict_spots(self):
        """Predict the spot locations on the image detector.
        
        The algorithm performs the following procedure:
        
         - First the set of miller indices are generated.
         
         - For each miller index, the rotation angle at which the diffraction
           conditions are met is calculated.
           
         - For those miller indices, i, where an angle cannot be calculated, 
           miller_indices[i] and rotation_angles[i] are removed from their 
           respective arrays
           
         - The rotation angles are then checked to see if they are within the
           rotation range.
           
         - For those reflections, i, not in the rotation range, miller_indces[i]
           and rotation_angles[i] are removed from their respective arrays.
           
         - The reciprocal lattice vectors are then calculated, followed by the
           diffracted beam vector for each reflection.
           
         - The image volume coordinates are then calculated for each reflection.
         
         - For each reflection, i, where the image volume coordinate couldnt be
           calculated, miller_indices[i], rotation_angles[i], beam_vectors[i] 
           and image_volume_coords[i] are removed from their respective arrays.
           
         - The image volume coordinates are then checked to see if they are
           within the image volume itself.
           
         - For those reflections, i, not in the image volume, miller_indices[i], 
           rotation_angles[i], beam_vectors[i] and image_volume_coords[i] are 
           removed from their respective arrays.
        
        """
        from dials.array_family import flex

        # Ensure internal arrays are reset
        self.reset()

        # Generate the miller indices
        self.__generate_miller_indices()
        
        # Calculate the set of valid rotation angles and miller indices and
        # remove the invalid miller indices and rotation angles and create 
        # an array of miller indices and rotation angles that correspond
        self.__rejig_miller_and_rotation_angle_arrays(
                self.__calculate_rotation_angles())

        # Filter the angles and miller indices with the rotation range
        self.__filter_reflections(self.__filter_angles_in_rotation_range())       
        
        # Calculate the beam vectors
        self._beam_vectors = self.__calculate_reciprocal_space_vectors() \
                           + self._beam.direction

        print self._goniometer

        from scitbx import matrix
        print (matrix.col(self._beam.direction).length(), 1.0 / self._beam.wavelength, 
            matrix.col(self._beam_vectors[0]).length())

        # Calculate the image volume coordinates and keep only those array 
        # elements that have a valid image coordinate
        self.__filter_reflections(self.__calculate_image_volume_coordinates())
 
        # Filter the image volume coordinates and remove any invalid spots to 
        # leave the valid ones remaining
        self.__filter_reflections(self.__filter_image_volume_coordinates())
         
    def __filter_array(self, array, keep):
        """Check the array is not null then remove the selected elements."""
        from dials.array_family import flex
        if array:
            return flex.remove_if_not(array, keep)
        else:
            return None

    def __filter_reflections(self, keep):
        """Remove elements from the reflecton arrays"""
        self._miller_indices  = self.__filter_array(self._miller_indices,  keep)
        self._rotation_angles = self.__filter_array(self._rotation_angles, keep)
        self._beam_vectors    = self.__filter_array(self._beam_vectors,    keep)
        self._image_coords    = self.__filter_array(self._image_coords,    keep)
  
    def __generate_miller_indices(self):
        """Generate a list of miller indices."""
        from dials.spot_prediction import IndexGenerator

        # Create the index generator
        index_generator = IndexGenerator(self._unit_cell, 
                                         self._space_group_type, 
                                         True, 
                                         self._d_min)
        
        # Generate and return the miller indices
        self._miller_indices = index_generator.to_array()

    def __calculate_rotation_angles(self):
        """Calculate the rotation angles for the given miller indices at which
        the diffracting condition is met.
        
        The function returns the rotation angles in an array of the form
        angles = flex.int(flex.grid(2, n)). The angles for a given miller index
        can be found at angles[0,i] and angles[1,i]

        Returns:
            
            Boolean array, True/False are angles valid
        
        """
        from dials.array_family import flex
        from dials.spot_prediction import RotationAngles
        from math import pi

        # Create the rotation angle calculator
        rotation_angle_calculator = RotationAngles(self._d_min, 
                                                   self._ub_matrix, 
                                                   self._beam.wavelength, 
                                                   self._goniometer.rotation_axis)

        # Calculate and return the rotation angles
        status = flex.bool()
        self._rotation_angles = rotation_angle_calculator.calculate(
            self.miller_indices, status)

        # Convert to degrees
        r2d = 180.0 / pi
        for i in range(len(self._rotation_angles)):
            self._rotation_angles[i] = self._rotation_angles[i] * r2d % 360

        # Return the status
        return status

    def __rejig_miller_and_rotation_angle_arrays(self, is_valid_angle):
        """Rejig the miller indices and rotation angles arrays to remove
        any elements corresponding to back rotation angles. Then rejig the
        arrays to make sure that for each element the miller index corresponds
        to the rotation angle.
        
        Args:
            is_valid_angle Is the angle valid (boolean array)
               
        """
        from dials.array_family import flex
        
        # Remove the invalid angles from the miller index array and add the 
        # array onto the end of itself.
        self._miller_indices = flex.remove_if_not(self.miller_indices, 
                                                  is_valid_angle)
        self._miller_indices = self.miller_indices.concatenate(self.miller_indices)

        # Remove rotation angles which are invalid and concatenate the arrays
        # so that we have a 1d array where for each element the angle matches
        # the miller indices at the same element in miller_indices
        rotation_angles_a = flex.remove_if_not(
            self.rotation_angles[0:1,:].as_1d(), is_valid_angle)
        rotation_angles_b = flex.remove_if_not(
            self.rotation_angles[1:2,:].as_1d(), is_valid_angle)
        self._rotation_angles = rotation_angles_a.concatenate(rotation_angles_b)

    def __filter_angles_in_rotation_range(self):
        """Check if each angle is within the rotation range.
                
        Returns:
            A boolean array, are the angles in the rotation range
        
        """
        from dials.spot_prediction import is_angle_in_range

        # Calculate the rotation angle range
        rotation_range = (self._goniometer.starting_angle,
                          self._goniometer.get_angle_from_frame(
                            self._goniometer.starting_frame + 
                            self._goniometer.num_frames))        

        # Check which angles are in the rotation range
        return is_angle_in_range(self.rotation_angles, rotation_range, deg=True)

    def __calculate_reciprocal_space_vectors(self):
        """Calculate the reciprocal space vectors of the given miller indices.
                    
        Returns:
            An array of reciprocal space vectors
        
        """
        from dials.geometry.transform import FromHklToRsv

        # Construct the hkl -> p transform 
        from_hkl_to_rsv = FromHklToRsv(self._ub_matrix, 
                                       self._goniometer.rotation_axis)
        
        # Calculate and return the set of reciprocal space vectors
        return from_hkl_to_rsv.apply(self.miller_indices, self.rotation_angles)

    def __calculate_image_volume_coordinates(self):
        """Calculate the image volume coordinates for the given beam vectors.
        
        The z coordinate is zero-based
                    
        Returns:
            An boolean array containing True/False status for each element
        
        """
        from dials.array_family import flex
        from dials.geometry.transform import FromBeamVectorToImageVolume

        # Construct the s1 -> xyz transform
        from_beam_vector_to_image_volume = FromBeamVectorToImageVolume(
                                            self._detector, self._goniometer)
        
        # Calculate and return the image volume coordinates
        status = flex.bool()
        self._image_coords = from_beam_vector_to_image_volume.apply(
                self.beam_vectors, self.rotation_angles, status)

        # Return the status
        return status

    def __filter_image_volume_coordinates(self):
        """Check if the image volume coordinates are within the image volume.
        
        Returns:
            A boolean array, are the coordinates within the image volume.
        
        """
        from dials.spot_prediction import in_volume

        # Calculate the x, y and z ranges
        range_x = (0, self._detector.size[0])
        range_y = (0, self._detector.size[1])
        range_z = (0, self._goniometer.num_frames)

        # Return an array in volume True/False
        return in_volume(self._image_coords, range_x, range_y, range_z)

