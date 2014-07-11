


#ifndef DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/integration/profile/grid_sampler_2d.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::Shoebox;

  class Flattened2DProfileFitting {
  public:

    typedef af::versa< double, af::c_grid<2> > double2d;
    typedef af::versa< int,    af::c_grid<2> > int2d;
    typedef af::ref< double, af::c_grid<2> > double2d_ref;
    typedef af::ref< int,    af::c_grid<2> > int2d_ref;

    Flattened2DProfileFitting(
        std::size_t image_width,
        std::size_t image_height,
        const af::const_ref< vec3<double> > &xyz,
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

      // Create the grid
      int2 image_size(image_width, image_height);
      int2 grid_size(5, 5);
      GridSampler2D grid(image_size, grid_size);

      // Interpolate shoeboxes to centre on pixel
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        interpolate(xyz[i], flattened_data[i].ref());
      }

      // Find the nearest grid point to each profile
      af::shared<std::size_t> profile_index(sbox.size());
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        double2 xy(xyz[i][0], xyz[i][1]);
        profile_index[i] = grid.nearest(xy);
      }

    }

    af::shared<double> intensity() const {
      return intensity_;
    }

    af::shared<double> variance() const {
      return variance_;
    }

  private:

    void interpolate(vec3<double> xyz, double2d_ref data) const {
      double x = xyz[0] - floor(xyz[0]) - 0.5;
      double y = xyz[1] - floor(xyz[1]) - 0.5;
      std::size_t h = data.accessor()[0];
      std::size_t w = data.accessor()[1];
      for (std::size_t j = 0; j < h; ++j) {
        for (std::size_t i = 0; i < w; ++i) {
          double f00 = data(j,i);
          double f10 = data(j+1,i);
          double f01 = data(j,i+1);
          double f11 = data(j+1,i+1);
          //I = f00*(1-x)*(1-y)+f10*x(1-y)+f01*(1-x)+f11*x*y;
        }
      }
    }

    af::shared<double> intensity_;
    af::shared<double> variance_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_FLATTENED_2D_PROFILE_FITTING_H
