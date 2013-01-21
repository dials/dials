
class SpotPredictor:

    def __init__(self, beam, goniometer, detector, n_image_frames, 
                 unit_cell, space_group_type, d_min, ub_matrix):

        # Set the equipment
        self.beam = beam
        self.goniometer = goniometer
        self.detector = detector
        self.n_image_frames = n_image_frames

        # Set the crystal stuff
        self.unit_cell = unit_cell
        self.space_group_type = space_group_type
        self.d_min = d_min
        self.ub_matrix = ub_matrix

        # Set all arrays to empty
        self.miller_indices = []
        self.rotation_angles = []
        self.beam_vectors = []
        self.image_volume_coords = []
    
    def predict_spots(self):

        from dials.spot_prediction import remove_if_not

        # Set all arrays to empty
        self.miller_indices = []
        self.rotation_angles = []
        self.beam_vectors = []
        self.image_volume_coords = []

        # Generate the miller indices
        miller_indices = self._generate_miller_indices()
        
        # Calculate the set of valid rotation angles and miller indices
        miller_indices, rotation_angles = self._calculate_rotation_angles(miller_indices)
        
        #reciprocal_space_vectors = self.unit_cell.reciprocal_space_vector(miller_indices)
        #print "1"
        #print self.ub_matrix
        #print self.unit_cell.reciprocal_parameters()
        #from scitbx.array_family import flex
        #from scitbx import matrix
        #beam_vectors = flex.vec3_double(len(reciprocal_space_vectors))
        #for i in range(0, len(reciprocal_space_vectors)):
        #    beam_vectors[i] = matrix.col(reciprocal_space_vectors[i]).rotate(matrix.col(self.goniometer.rotation_axis), rotation_angles[i]) + matrix.col(self.beam.direction)
        #from cctbx.uctbx import unit_cell
        #print "a"
        #print self.ub_matrix

        #print self.unit_cell.fractionalization_matrix()
        
        from dials.geometry.transform import FromHklToRsv

        from_hkl_to_rsv = FromHklToRsv(self.ub_matrix, self.goniometer.rotation_axis)
        rsv = from_hkl_to_rsv.apply(miller_indices, rotation_angles)

        from scitbx.array_family import flex
        beam_vectors = rsv + self.beam.direction

        # Calculate the beam vectors
        #beam_vectors = self._calculate_beam_vectors(miller_indices, rotation_angles)

        # Calculate the image volume coordinates
        status = flex.bool(0)
        image_volume_coords = self._calculate_image_volume_coordinates(beam_vectors, rotation_angles, status)

        self.miller_indices      = remove_if_not(miller_indices,      status)
        self.rotation_angles     = remove_if_not(rotation_angles,     status)
        self.beam_vectors        = remove_if_not(beam_vectors,        status)
        self.image_volume_coords = remove_if_not(image_volume_coords, status)

        # Filter the image volume coordinates
        is_valid_coord = self._filter_image_volume_coordinates(image_volume_coords)

        # Remove any invalid spots to leave the valid ones remaining
        self.miller_indices      = remove_if_not(miller_indices,      is_valid_coord)
        self.rotation_angles     = remove_if_not(rotation_angles,     is_valid_coord)
        self.beam_vectors        = remove_if_not(beam_vectors,        is_valid_coord)
        self.image_volume_coords = remove_if_not(image_volume_coords, is_valid_coord)
 
    def _generate_miller_indices(self):
        from dials.spot_prediction import IndexGenerator

        # Create the index generator
        index_generator = IndexGenerator(self.unit_cell, 
                                         self.space_group_type, 
                                         True, 
                                         self.d_min)
        
        # Generate and return the miller indices
        return index_generator.to_array()

    def _calculate_rotation_angles(self, miller_indices):
         
        from dials.spot_prediction import RotationAngles, is_angle_in_range
        from dials.spot_prediction import remove_if_not

        # Calculate the rotation angle range
        rotation_angle_range = (self.goniometer.starting_angle,
                                self.goniometer.get_angle_from_frame(
                                    self.n_image_frames + 
                                    self.goniometer.starting_frame))

        # Create the rotation angle calculator
        rotation_angle_calculator = RotationAngles(self.d_min, 
                                                   self.ub_matrix, 
                                                   self.beam.wavelength, 
                                                   self.goniometer.rotation_axis)

        # Calculate the rotation angles
        from scitbx.array_family import flex
        status = flex.bool(0)
        rotation_angles = rotation_angle_calculator.calculate(miller_indices, status)
 
        miller_indices = remove_if_not(miller_indices, status)
        rotation_angles_a = remove_if_not(rotation_angles[0:1,:].as_1d(), status)
        rotation_angles_b = remove_if_not(rotation_angles[1:2,:].as_1d(), status)
  
        miller_indices = miller_indices.concatenate(miller_indices)
        rotation_angles = rotation_angles_a.concatenate(rotation_angles_b)
  
        # Get the corresponding miller indices
        #miller_indices = rotation_angle_calculator.miller_indices()


        in_rotation_range = is_angle_in_range(rotation_angles, rotation_angle_range)
        
        rotation_angles = remove_if_not(rotation_angles, in_rotation_range)
        miller_indices = remove_if_not(miller_indices, in_rotation_range)

        # Return the rotation angles and miller indices as a tuple
        return (miller_indices, rotation_angles)

    def _calculate_beam_vectors(self, miller_indices, rotation_angles):

        from dials.geometry.transform import FromHklToBeamVector

        # Construct the hkl -> s1 transform 
        from_hkl_to_beam_vector = FromHklToBeamVector(
                                    self.ub_matrix, 
                                    self.beam.direction, 
                                    self.goniometer.rotation_axis)
        
        # Calculate and return the set of beam vectors
        return from_hkl_to_beam_vector.apply(miller_indices, rotation_angles)

    def _calculate_image_volume_coordinates(self, beam_vectors, rotation_angles, status):

        from dials.geometry.transform import FromBeamVectorToImageVolume

        # Construct the s1 -> xyz transform
        from_beam_vector_to_image_volume = FromBeamVectorToImageVolume(
                                                self.detector, 
                                                self.goniometer)
        
        # Calculate and return the image volume coordinates
        return from_beam_vector_to_image_volume.apply(beam_vectors, 
                                                      rotation_angles,
                                                      status)

    def _filter_image_volume_coordinates(self, image_volume_coords):
            
        from dials.spot_prediction import in_volume

        # Calculate the x, y and z ranges
        range_x = (0, self.detector.size[0])
        range_y = (0, self.detector.size[1])
        range_z = (0, self.n_image_frames)

        # Return an array in volume True/False
        return in_volume(image_volume_coords, range_x, range_y, range_z)
