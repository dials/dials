
class IndexGenerator:
    """A class to generate hkl indices from a unit cell."""

    def __init__(self, uc, cell_space_group, dmin):
        """Generate all the possible reflection indices from a unit cell object
        and calculate the symmetry reduced set of reflections. The list of
        reflection indices are then saved in self.indices.

        Args:
            uc: The unit cell object
            dmin: The resolution
            cell_space_group: The space group object

        """
        indices = self._generate_reflection_indices(uc, dmin)
        #self.indices = indices
        self.indices = self._remove_absent_indices(indices, cell_space_group)

    def _generate_reflection_indices(self, uc, dmin):
        """Generate the possible reflection indices from a unit cell object:
        N.B. that these are not symmetry reduced.

        Args:
            uc: The unit cell object
            dmin: The resolution

        Returns:
            The List of indices

        """
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

    def _remove_absent_indices(self, indices, cell_space_group):
        """From the given list of indices, remove those reflections which
        should be systematic absences according to the given space group.

        Args:
            indices: The generated indices
            cell_space_group: The space group object

        Returns:
            The symmetry reduced indices

        """
        # The list of present reflections
        present = []

        # Check that each reflection is present
        for hkl in indices:
            if not cell_space_group.is_sys_absent(hkl):
                present.append(hkl)

        # Return those reflections that are present
        return present
