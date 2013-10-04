/*
 *  mask_builder.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MASK_BUILDER_H
#define DIALS_MASK_BUILDER_H
#include <iostream>
#include <cmath>
#include <dials/array_family/scitbx_shared_and_versa.h>


namespace dials { namespace algorithms { namespace shoebox {

  af::versa< int, af::c_grid<2> > build_mask(int nx, int ny, int nrx,
                                                 int nry, int nc) {

    std::cout << "\n" << "nx=" << nx <<"\n" << "ny=" << ny <<"\n"
    << "nrx=" << nrx <<"\n" << "nry=" << nry <<"\n" << "nc=" << nc <<"\n";

    // creating a versa array the size of the entered parameters
    // with all pixels assigned to valid but not in peak area (= 3)
    af::versa< int, af::c_grid<2> > mask(af::c_grid<2>(ny, nx),3);



    /*

            <------------------ NX = 23 --------------->

       ^    - - - - - - - - - - - - - - - - - - - - - - -  ^
       !    - - - - - - - - - - - - - - - - - - - - - - -    NRY =2
       !    - - - - - -                       - - - - - -  ^
       !    - - - - -                           - - - - -
       !    - - - -                               - - - -
       !    - - -                                   - - -
       !    - - -                                   \ - -
       !    - - -                                   - \ -
    NY  =17 - - -                                   - - \___
       !    - - -                                   - - /   ^
       !    - - -                                   - / -   !
       !    - - -                                   / - -   !
       !    - - - -                               - - - -   !
       !    - - - - -                           - - - - -   NC =8
       !    - - - - - -                       - - - - - -   !
       !    - - - - - - - - - - - - - - - - - - - - - - -   !
       ^    - - - - - - - - - - - - - - - - - - - - - - -   ^
            <NRX> = 3

    */

    // looping through pixels and filtering conditions saved
    // in the parameters nc, nrx and nry
    for(int row = 0; row < ny; row++){
      for( int col = 0; col < nx; col++ ){

        // conditioning first the corners an then the boundaries
          if( row + col >= nc and row < col + ny - nc and
           row > col - nx + nc and (ny-row)+(nx-col) > nc + 1  and
            row >= nry and row  < (ny - nry) and
             col >= nrx and col  < (nx - nrx) )
          {
              // when [ pixel mask ] = 5 both bits are set:
              // is valid and is in peak area
              mask(row,col) = 5;
          }
      }
    }
    return mask;

  }
}}}

#endif
