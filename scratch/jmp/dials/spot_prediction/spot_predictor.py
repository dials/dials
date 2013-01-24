
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
        return self._image_volume_coords
    
    def reset(self):
        """Reset the internal arrays to None."""
        self._miller_indices = None
        self._rotation_angles = None
        self._beam_vectors = None
        self._image_volume_coords = None
    
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
        miller_indices = self.__generate_miller_indices()
        
        # Calculate the set of valid rotation angles and miller indices and
        # remove the invalid miller indices and rotation angles and create 
        # an array of miller indices and rotation angles that correspond
        is_valid_angle = flex.bool()
        rotation_angles = self.__calculate_rotation_angles(miller_indices, is_valid_angle)
        miller_indices = flex.remove_if_not(miller_indices, is_valid_angle)
        rotation_angles_a = flex.remove_if_not(rotation_angles[0:1,:].as_1d(), is_valid_angle)
        rotation_angles_b = flex.remove_if_not(rotation_angles[1:2,:].as_1d(), is_valid_angle)
        miller_indices = miller_indices.concatenate(miller_indices)
        rotation_angles = rotation_angles_a.concatenate(rotation_angles_b)        
        
        # Filter the angles and miller indices with the rotation range
        in_rotation_range = self.__filter_angles_in_rotation_range(rotation_angles)
        rotation_angles = flex.remove_if_not(rotation_angles, in_rotation_range)
        miller_indices  = flex.remove_if_not(miller_indices, in_rotation_range)        
        
        # Calculate the beam vectors
        beam_vectors = self.__calculate_reciprocal_space_vectors(
                            miller_indices, rotation_angles) + self._beam.direction

        # Calculate the image volume coordinates and keep only those array 
        # elements that have a valid image coordinate
        is_valid_coord = flex.bool()
        image_volume_coords = self.__calculate_image_volume_coordinates(
                                beam_vectors, rotation_angles, is_valid_coord)
        miller_indices = flex.remove_if_not(miller_indices, is_valid_coord)
        rotation_angles = flex.remove_if_not(rotation_angles, is_valid_coord)
        beam_vectors = flex.remove_if_not(beam_vectors, is_valid_coord)
        image_volume_coords = flex.remove_if_not(image_volume_coords, is_valid_coord)

        # Filter the image volume coordinates and remove any invalid spots to 
        # leave the valid ones remaining
        is_valid_coord = self.__filter_image_volume_coordinates(image_volume_coords)
        self._miller_indices = flex.remove_if_not(miller_indices, is_valid_coord)
        self._rotation_angles = flex.remove_if_not(rotation_angles, is_valid_coord)
        self._beam_vectors = flex.remove_if_not(beam_vectors, is_valid_coord)
        self._image_volume_coords = flex.remove_if_not(image_volume_coords, is_valid_coord)
 
    def __generate_miller_indices(self):
        """Generate a list of miller indices.
        
        Returns:
            An array of miller indices
        
        """
        from dials.spot_prediction import IndexGenerator

        # Create the index generator
        index_generator = IndexGenerator(self._unit_cell, 
                                         self._space_group_type, 
                                         True, 
                                         self._d_min)
        
        # Generate and return the miller indices
        return index_generator.to_array()

    def __calculate_rotation_angles(self, miller_indices, status):
        """Calculate the rotation angles for the given miller indices at which
        the diffracting condition is met.
        
        The function returns the rotation angles in an array of the form
        angles = flex.int(flex.grid(2, n)). The angles for a given miller index
        can be found at angles[0,i] and angles[1,i]
        
        Args:
            miller_indices The array of miller indices
            status Boolean array, were rotation angles calculated
            
        Returns:
            
            A 2xn array of rotation angles .
        
        """
        from dials.spot_prediction import RotationAngles

        # Create the rotation angle calculator
        rotation_angle_calculator = RotationAngles(self._d_min, 
                                                   self._ub_matrix, 
                                                   self._beam.wavelength, 
                                                   self._goniometer.rotation_axis)

        # Calculate and return the rotation angles
        return rotation_angle_calculator.calculate(miller_indices, status)

    def __filter_angles_in_rotation_range(self, rotation_angles):
        """Check if each angle is within the rotation range.
        
        Args:
            rotation_angles The list of rotation angles
        
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
        return is_angle_in_range(rotation_angles, rotation_range)

    def __calculate_reciprocal_space_vectors(self, miller_indices, 
                                             rotation_angles):
        """Calculate the reciprocal space vectors of the given miller indices.
        
        Args:
            miller_indices The array of miller indices
            rotation_angles The array of rotation angles
            
        Returns:
            An array of reciprocal space vectors
        
        """
        from dials.geometry.transform import FromHklToRsv

        # Construct the hkl -> p transform 
        from_hkl_to_rsv = FromHklToRsv(self._ub_matrix, 
                                       self._goniometer.rotation_axis)
        
        # Calculate and return the set of reciprocal space vectors
        return from_hkl_to_rsv.apply(miller_indices, rotation_angles)

    def __calculate_image_volume_coordinates(self, beam_vectors, 
                                             rotation_angles, status):
        """Calculate the image volume coordinates for the given beam vectors.
        
        The z coordinate is zero-based
        
        Args:
            beam_vectors The array of beam vectors
            rotation_angles The array of rotation angles
            status Boolean array, was calculation successful
            
        Returns:
            An array of (x, y, z) image volume coordinates
        
        """
        from dials.geometry.transform import FromBeamVectorToImageVolume

        # Construct the s1 -> xyz transform
        from_beam_vector_to_image_volume = FromBeamVectorToImageVolume(
                                            self._detector, self._goniometer)
        
        # Calculate and return the image volume coordinates
        return from_beam_vector_to_image_volume.apply(
                beam_vectors, rotation_angles, status)

    def __filter_image_volume_coordinates(self, image_volume_coords):
        """Check if the image volume coordinates are within the image volume.
                
        Args:
            image_volume_coords The image volume coordinates
        
        Returns:
            A boolean array, are the coordinates within the image volume.
        
        """
        from dials.spot_prediction import in_volume

        # Calculate the x, y and z ranges
        range_x = (0, self._detector.size[0])
        range_y = (0, self._detector.size[1])
        range_z = (0, self._goniometer.num_frames)

        # Return an array in volume True/False
        return in_volume(image_volume_coords, range_x, range_y, range_z)

