class centroid_interface_prototype( object ):
    '''Centroid calculation interface: all of the stuff in here is in pixels.'''

    def __init__( self, bounding_boxes, dxtbx_sweep_object ):
        self._bounding_boxes = bounding_boxes
        self._sweep = dxtbx_sweep_object
        self._compute_centroids()

        return

    def get_centroids( self ):
        return self._centroids

    def compute_centroid_from_bbox( self, bbox ):
        '''Overload me: bbox has form tuple:

        (fmin, fmax, rmin, rmax, cmin, cmax)

        where f is frame #, r is row #, c is column #. Should return f, r, c
        in pixels of the centroid position, and df, dr, dc standard deviations
        of the same.'''

        raise RuntimeError, 'overload me'

    def _compute_centroids( self ):

        self._centroids = { }

        for hkl in self._bounding_boxes:

            self._centroids[hkl] = []

            for bbox in self._bounding_boxes[hkl]:
                self._centroids[hkl].append( self.compute_centroid_from_bbox( 
                    bbox ) )

        return



