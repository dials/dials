


#ifndef DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::Shoebox;

  class Flattened2DProfileFitting {
  public:

    typedef af::versa< double, af::c_grid<2> > double2d;
    typedef af::versa< int,    af::c_grid<2> > int2d;

    Flattened2DProfileFitting(
        const af::const_ref< Shoebox<> > &sbox) {

      // Allocate arrays for flattened shoebox
      af::shared<double2d> flattened_data(sbox.size());
      af::shared<double2d> flattened_bgrd(sbox.size());
      af::shared<int2d>    flattened_mask(sbox.size());

      // Create the flattened shoeboxes
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        af::c_grid<2> grid(sbox[i].ysize(), sbox[i].xsize());
        flattened_data[i] = double2d(grid, 0);
        flattened_bgrd[i] = double2d(grid, 0);
        flattened_mask[i] = int2d(grid, 0);
        for (std::size_t z = 0; z < sbox[i].zsize(); ++z) {
          for (std::size_t y = 0; y < sbox[i].ysize(); ++y) {
            for (std::size_t x = 0; x < sbox[i].xsize(); ++x) {
              flattened_data[i](y,x) += sbox[i].data(z,y,x);
              flattened_bgrd[i](y,x) += sbox[i].background(z,y,x);
              flattened_mask[i](y,x) |= sbox[i].mask(z,y,x);
            }
          }
        }
      }
    }

    af::shared<double> intensity() const {
      return intensity_;
    }

    af::shared<double> variance() const {
      return variance_;
    }

  private:

    af::shared<double> intensity_;
    af::shared<double> variance_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
