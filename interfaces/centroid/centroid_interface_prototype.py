class centroid_interface_prototype( object ):
    '''Centroid calculation interface: all of the stuff in here is in pixels.'''

    def __init__(self, reflections):
        self._reflections = reflections
        self._compute_centroids()

        return

    def get_reflections( self ):
        return self._reflections

    def compute_centroid_from_bbox( self, image ):
        '''Overload me: bbox has form tuple:

        (fmin, fmax, rmin, rmax, cmin, cmax)

        where f is frame #, r is row #, c is column #. Should return f, r, c
        in pixels of the centroid position, and df, dr, dc standard deviations
        of the same.'''

        raise RuntimeError, 'overload me'

    def _compute_centroids( self ):

        self._centroids = { }

        for i, ref in enumerate(self._reflections):

            # Calculate the centroid
            f, r, c, sf, sr, sc = self.compute_centroid_from_bbox(ref.image)

            # Add shoebox offset to centroid
            shoebox = ref.shoebox
            f += ref.shoebox[4]
            r += ref.shoebox[2]
            c += ref.shoebox[0]

            # Add centroid data to reflection
            ref.centroid_position = (c, r, f)
            ref.centroid_variance = (sc, sr, sf)

            # Copy reflection back into array
            self._reflections[i] = ref

        return
