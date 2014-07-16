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

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/integration/profile/grid_sampler_2d.h>
#include <dials/algorithms/integration/interpolate_profile2d.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::Shoebox;

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



    }

    af::shared<std::size_t> index() {
      return index_;
    }

    af::shared<int6> bbox() {
      return bbox_;
    }

    af::shared<double2d> data() {
      return data_;
    }

    af::shared<double2d> background() {
      return bgrd_;
    }

    af::shared<int2d> mask() {
      return mask_;
    }

    af::shared<std::size_t> index_;
    af::shared<double2d> data_;
    af::shared<double2d> bgrd_;
    af::shared<int2d> mask_;
    af::shared<int6> bbox_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_FLATTENER_H
