//+++ b/trunk/scratch/luiso_s/lui_example_helper.h

#ifndef DIALS_MASK_BUILDER_H
#define DIALS_MASK_BUILDER_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>
#include <cmath>
//#include <dials/array_family/scitbx_shared_and_versa.h>


namespace dials { namespace algorithms { namespace shoebox {
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  //af::versa< double, af::c_grid<2> > smooth_2d(
  flex_int build_mask(int nx, int ny, int nrx, int nry, int nc) {

    std::cout << "\n" << "nx=" << nx <<"\n"
              << "ny="  << ny  <<"\n"
              << "nrx=" << nrx <<"\n"
              << "nry=" << nry <<"\n"
              << "nc="  << nc  <<"\n";
    //af::versa< double, af::c_grid<2> > data2dsmoth(data2d.accessor(),0);
    flex_int mask(flex_grid<>(ny, nx),3);

    return mask;

  }
}}}

#endif
