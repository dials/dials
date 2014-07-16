/*
 * shoebox_flattener.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_FLATTENER_H
#define DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_FLATTENER_H

#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/integration/profile/grid_sampler_2d.h>
#include <dials/algorithms/integration/interpolate_profile2d.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::tiny;
  using model::Shoebox;

  /**
   * Crude class to flatten and resize shoeboxes in each grid area
   */
  class ShoeboxFlattener {
  public:

    typedef af::versa< double, af::c_grid<2> > double2d;
    typedef af::versa< int,    af::c_grid<2> > int2d;

    ShoeboxFlattener(
        const GridSampler2D &grid,
        const af::const_ref< vec3<double> > &xyz,
        const af::const_ref< Shoebox<> > &sbox)
          : index_(xyz.size()),
            data_(xyz.size()),
            bgrd_(xyz.size()),
            mask_(xyz.size()),
            bbox_(xyz.size()) {
      DIALS_ASSERT(xyz.size() == sbox.size());

      // Get the list of grid indices
      for (std::size_t i = 0; i < index_.size(); ++i) {
        index_[i] = grid.nearest(double2(xyz[i][0], xyz[i][1]));
      }

      // Find the minimum sizes of the shoeboxes for each grid
      af::shared<std::size_t> count(xyz.size(), 0);
      af::shared< tiny<std::size_t,4> > min_size(
          xyz.size(), tiny<std::size_t,4>(0, 0, 0, 0));
      for (std::size_t i = 0; i < index_.size(); ++i) {
        std::size_t j = index_[i];
        int x0 = sbox[i].bbox[0];
        int x1 = sbox[i].bbox[1];
        int y0 = sbox[i].bbox[2];
        int y1 = sbox[i].bbox[3];
        int x = (int)floor(xyz[i][0]);
        int y = (int)floor(xyz[i][1]);
        DIALS_ASSERT(x1 > x0 && y1 > y0);
        DIALS_ASSERT(x > x0 && x1 > x);
        DIALS_ASSERT(y > y0 && y1 > y);
        std::size_t xsize1 = x - x0;
        std::size_t xsize2 = x1 - x;
        std::size_t ysize1 = y - y0;
        std::size_t ysize2 = y1 - y;
        if (count[j] == 0) {
          min_size[j][0] = xsize1;
          min_size[j][1] = xsize2;
          min_size[j][2] = ysize1;
          min_size[j][3] = ysize2;
        } else {
          if (xsize1 < min_size[j][0]) min_size[j][0] = xsize1;
          if (xsize2 < min_size[j][1]) min_size[j][1] = xsize2;
          if (ysize1 < min_size[j][2]) min_size[j][2] = ysize1;
          if (ysize2 < min_size[j][3]) min_size[j][3] = ysize2;
        }
        count[j]++;
      }

      // Create flattened shoeboxes
      for (std::size_t i = 0; i < index_.size(); ++i) {
        std::size_t j = index_[i];
        std::size_t xsize1 = min_size[j][0];
        std::size_t xsize2 = min_size[j][1];
        std::size_t ysize1 = min_size[j][2];
        std::size_t ysize2 = min_size[j][3];
        DIALS_ASSERT(xsize1 > 0);
        DIALS_ASSERT(xsize2 > 0);
        DIALS_ASSERT(ysize1 > 0);
        DIALS_ASSERT(ysize2 > 0);
        int x = (int)floor(xyz[i][0]);
        int y = (int)floor(xyz[i][1]);
        int x0 = x - xsize1;
        int x1 = x + xsize2;
        int y0 = y - ysize1;
        int y1 = y + ysize2;
        bbox_[i][0] = x0;
        bbox_[i][1] = x1;
        bbox_[i][2] = y0;
        bbox_[i][3] = y1;
        bbox_[i][4] = sbox[i].bbox[4];
        bbox_[i][5] = sbox[i].bbox[5];
        std::size_t xsize = xsize1 + xsize2;
        std::size_t ysize = ysize1 + ysize2;
        int x00 = sbox[i].bbox[0];
        int x11 = sbox[i].bbox[1];
        int y00 = sbox[i].bbox[2];
        int y11 = sbox[i].bbox[3];
        int z00 = sbox[i].bbox[4];
        int z11 = sbox[i].bbox[5];
        DIALS_ASSERT(x11 > x00);
        DIALS_ASSERT(y11 > y00);
        DIALS_ASSERT(z11 > z00);
        std::size_t old_xsize = x11 - x00;
        std::size_t old_ysize = y11 - y00;
        std::size_t old_zsize = z11 - z00;
        DIALS_ASSERT(old_xsize >= xsize);
        DIALS_ASSERT(old_ysize >= ysize);
        DIALS_ASSERT(x0 >= x00);
        DIALS_ASSERT(y0 >= y00);
        std::size_t xoff = x0 - x00;
        std::size_t yoff = y0 - y00;
        double2d data(af::c_grid<2>(ysize, xsize),0);
        double2d bgrd(af::c_grid<2>(ysize, xsize),0);
        int2d mask(af::c_grid<2>(ysize, xsize),0);
        for (std::size_t kk = 0; kk < old_zsize; ++kk) {
          for (std::size_t jj = 0; jj < ysize; ++jj) {
            for (std::size_t ii = 0; ii < xsize; ++ii) {
              data(jj,ii) += sbox[i].data(kk,jj+yoff,ii+xoff);
              bgrd(jj,ii) += sbox[i].background(kk,jj+yoff,ii+xoff);
              mask(jj,ii) |= sbox[i].mask(kk,jj+yoff,ii+xoff);
            }
          }
        }
        data_[i] = data;
        bgrd_[i] = bgrd;
        mask_[i] = mask;
      }
    }

    af::shared<std::size_t> index() {
      return index_;
    }

    af::shared<int6> bbox() {
      return bbox_;
    }

    double2d data(std::size_t index) {
      return data_[index];
    }

    double2d background(std::size_t index) {
      return bgrd_[index];
    }

    int2d mask(std::size_t index) {
      return mask_[index];
    }

    std::size_t size() const {
      return index_.size();
    }

    af::shared<std::size_t> index_;
    af::shared<double2d> data_;
    af::shared<double2d> bgrd_;
    af::shared<int2d> mask_;
    af::shared<int6> bbox_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_FLATTENER_H
