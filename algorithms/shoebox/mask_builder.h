//+++ b/trunk/scratch/luiso_s/lui_example_helper.h

#ifndef DIALS_MASK_BUILDER_H
#define DIALS_MASK_BUILDER_H
#include <iostream>
#include <cmath>
#include <dials/array_family/scitbx_shared_and_versa.h>


namespace dials { namespace algorithms { namespace shoebox {

  af::versa< int, af::c_grid<2> > build_mask(int nx, int ny, int nrx,
                                                 int nry, int nc) {

    std::cout << "\n" << "nx=" << nx <<"\n"
                      << "ny=" << ny <<"\n"
                    << "nrx=" << nrx <<"\n"
                    << "nry=" << nry <<"\n"
                    << "nc="  << nc  <<"\n";

    af::versa< int, af::c_grid<2> > mask(af::c_grid<2>(ny, nx),3);
    for(int row = 0; row < ny; row++){
      for( int col = 0; col < nx; col++ ){
          if( row + col >= nc and row < col + ny - nc and
           row > col - nx + nc and (ny-row)+(nx-col) > nc + 1  and
            row >= nry and row  < (ny - nry) and
            col >= nrx and col  < (nx - nrx) )
          {
          mask(row,col) = 5;
        }
      }
    }
    return mask;

  }
}}}

#endif
