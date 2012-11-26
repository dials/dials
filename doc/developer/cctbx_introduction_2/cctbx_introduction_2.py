#!/usr/bin/env cctbx.python
# cctbx_introduction_2.py
#
# An illustration of how to use cctbx code with orientation matrices, and
# how to perform a simple geometric strategy calculation, thus demonstrating
# summetry...

from scitbx import matrix
from cctbx import uctbx

def parse_mosflm_matrix(matrix_file):
    '''Parse the mosflm matrix file to get: U, B, wavelength, unit cell.'''

    tokens = open(matrix_file, 'r').read(512).replace('-', ' -').split()

    assert(len(tokens) in (30, 32))

    lUB = matrix.sqr(map(float, tokens[0:9]))
    U = matrix.sqr(map(float, tokens[12:21]))
    
    uc = uctbx.unit_cell(map(float, tokens[21:27]))

    # derive the wavelength

    A = lUB.inverse().elems
    a = matrix.col(A[0:3])
    b = matrix.col(A[3:6])
    c = matrix.col(A[6:9])

    wavelength = (uc.parameters()[0] / a.length() +
                  uc.parameters()[1] / b.length() +
                  uc.parameters()[2] / c.length()) / 3

    # now verify that we reproduce the lUB matrix

    B = matrix.sqr(uc.fractionalization_matrix())

    m_format = '%7.4f %7.4f %7.4f\n%7.4f %7.4f %7.4f\n%7.4f %7.4f %7.4f' 

    print m_format % B.elems
    print m_format % ((1.0 / wavelength) * U.inverse() * lUB).elems

if __name__ == '__main__':

    import sys
    
    parse_mosflm_matrix(sys.argv[1])
    

    
    
