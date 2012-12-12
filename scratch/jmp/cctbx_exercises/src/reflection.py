

class Generator:
    
    def __init__(self):
        pass
    
    def generate_reflection_indices(self, uc, dmin):
        """Generate the possible reflection indices from a unit cell object: N.B.
        that these are not symmetry reduced."""
       
        # Get the maximum hkl indices from the unit cell parameters
        maxh, maxk, maxl = uc.max_miller_indices(dmin)
    
        # The list of reflection indices
        indices = []
    
        # Loop through the maximum range of hkl indices
        for h in range(-maxh, maxh + 1):
            for k in range(-maxk, maxk + 1):
                for l in range(-maxl, maxl + 1):
                    
                    # Do not include (0, 0, 0)
                    if h == 0 and k == 0 and l == 0:
                        continue
                    
                    # Ensure the resolution is less than the minimum resolution
                    if uc.d((h, k, l)) < dmin:
                        continue
                    
                    # Append indices to list
                    indices.append((h, k, l))
    
        # Return the list of indices
        return indices


    def remove_absent_indices(self, indices, cell_space_group):
        """From the given list of indices, remove those reflections which should be
        systematic absences according to the given space group."""
       
        # The list of present reflections
        present = []
        
        # Check that each reflection is present
        for hkl in indices:
            if not cell_space_group.is_sys_absent(hkl):
                present.append(hkl)
    
        # Return those reflections that are present
        return present


    def generate_intersection_angles(self, a_matrix, dmin, wavelength, indices):
        """From an A matrix following the Mosglm convention and the list of indices,
        return a list of phi (hkl) where (typically) there will be 2 records
        corresponding to each h, k, l."""
        
        from rstbx.diffraction import rotation_angles
        from scitbx import matrix
        import math
        
        # Construct an object to calculate the rotation of a reflection about
        # the (0, 1, 0) axis.
        ra = rotation_angles(dmin, a_matrix, wavelength, matrix.col((1, 0, 0)))
            
        # Convert radians to degrees
        r2d = 180.0 / math.pi
    
        # Loop through all the given reflection indices    
        phi_hkl = []
        for hkl in indices:
            if ra(hkl):
                
                # Calculate the 2 intersection angles for the reflection index
                # and append them to the list of phi-angles
                phis = ra.get_intersection_angles()
                phi_hkl.append((phis[0] * r2d % 360, hkl))
                phi_hkl.append((phis[1] * r2d % 360, hkl))
                
        # Return the list of phi-angles
        return phi_hkl


    def generate_observed_reflection_indices(self, ub_matrix, unit_cell, cell_space_group, 
                                             dmin, wavelength):
        """Predict the reflections.
        
        Calculate the indices of all predicted reflections based on the unit cell 
        parameters and the given resolution. Then remove the indices that are
        absent due to the symmetry. Lastly, generate the intersection angles, phi,
        and return these.
        
        :param cbf_handle: The CBF file handle
        :param ub_matrix: The UB matrix
        :param unit_cell: The unit cell parameters
        :param cell_space_group: The unit cell space group symmetry
        :param dmin: The resolution
        :param wavelength: The beam wavelength
        :returns: A list of reflection indices
        
        """
        # Generate reflection indices from the unit cell parameters and resolution.
        # Then remove indices that are systemmatically absent because of symmetry.
        # Generate the intersection angles of the remaining reflections and return.
        indices = self.generate_reflection_indices(unit_cell, dmin)
        present = self.remove_absent_indices(indices, cell_space_group)
        return self.generate_intersection_angles(ub_matrix, dmin, wavelength, present)
    
    
    def _get_intersection_angles(self, h, ra):
        from math import pi
        r2d = 180.0 / pi
        if ra(h):
            return tuple(ra.get_intersection_angles()) * r2d % 360
        else:
            return tuple(0, 0)
    
    def generate(self, uc, dmin, a_matrix, wavelength):
        
        from rstbx.diffraction import rotation_angles
        from scitbx import matrix
        import math
        
        indices = uc.complete_miller_set_with_lattice_symmetry(False, dmin)

        # Construct an object to calculate the rotation of a reflection about
        # the (0, 1, 0) axis.
        ra = rotation_angles(dmin, a_matrix, wavelength, matrix.col((1, 0, 0)))
        angles = map(lambda h: self._get_intersection_angles(h, ra), indices)
    
        return (indices, angles)
    
if __name__ == '__main__':
    
    
    gen = Generator()
    