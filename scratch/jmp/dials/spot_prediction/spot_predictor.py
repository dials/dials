
class SpotPredictor:

    def __init__(self, beam, goniometer, detector, n_image_frames, 
                 unit_cell, space_group_type, d_min):

        # Set the equipment
        self.beam = beam
        self.goniometer = goniometer
        self.detector = detector
        self.n_image_frames = n_image_frames

        # Set the crystal stuff
        self.unit_cell = unit_cell
        self.space_group_type = space_group_type
        self.d_min

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

        # Calculate the beam vectors
        beam_vectors = self._calculate_beam_vectors(miller_indices, rotation_angles)

        # Calculate the image volume coordinates
        image_volume_coords = self._calculate_image_volume_coordinates(beam_vectors, rotation_angles)

        # Filter the image volume coordinates
        is_valid_coord = self._filer_image_volume_coordinates(image_volume_coords)

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
         
        from dials.spot_prediction import RotationAngles

        # Calculate the rotation angle range
        rotation_angle_range = (self.goniometer.starting_angle,
                                self.goniometer.from_frame_to_angle(
                                    self.n_image_frames + 
                                    self.goniometer.starting_frame))

        # Create the rotation angle calculator
        rotation_angle_calculator = RotationAngles(self.d_min, 
                                                   self.unit_cell.orthogonalization_matrix, 
                                                   self.beam.wavelength, 
                                                   self.beam.rotation_axis, 
                                                   rotation_angle_range)

        # Calculate the rotation angles
        rotation_angles = rotation_angle_calculator.calculate(miller_indices)
  
        # Get the corresponding miller indices
        miller_indices = rotation_angle_calculator.miller_indices()

        # Return the rotation angles and miller indices as a tuple
        return (miller_indices, rotation_angles)

    def _calculate_beam_vectors(self, miller_indices, rotation_angles):

        from dials.geometry import ReciprocalLatticeCoordinateSystem
        from dials.geometry.transform import FromHklToBeamVector

        # Create the reciprocal lattice coordinate system
        reciprocal_lattice_coordinate_system = ReciprocalLatticeCoordinateSystem(
                                                unit_cell.orthogonalization_matrix)

        # Construct the hkl -> s1 transform 
        from_hkl_to_beam_vector = FromHklToBeamVector(
                                    reciprocal_lattice_coordinate_system, 
                                    self.beam.direction, 
                                    self.goniometer.rotation_axis)
        
        # Calculate and return the set of beam vectors
        return from_hkl_to_beam_vector.apply(miller_indices, rotation_angles)

    def _calculate_image_volume_coordinates(self, beam_vectors, rotation_angles):

        from dials.geometry.transform import FromBeamVectorToImageVolume

        # Construct the s1 -> xyz transform
        from_beam_vector_to_image_volume = FromBeamVectorToImageVolume(
                                                self.detector, 
                                                self.goniometer)
        
        # Calculate and return the image volume coordinates
        return from_beam_vector_to_image_volume.apply(beam_vectors, 
                                                      rotation_angles)

    def _filter_image_volume_coordinates(self, image_volume_coords):
            
        from dials.spot_prediction import in_volume

        # Calculate the x, y and z ranges
        range_x = (0, self.detector.size[0] - 1)
        range_y = (0, self.detector.size[1] - 1)
        range_z = (0, self.n_image_frames - 1)

        # Return an array in volume True/False
        return in_volume(image_volume_coords, range_x, range_y, range_z)