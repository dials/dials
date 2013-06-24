/*
 * grid_sampler.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_H

#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::double2;

  /**
   * Class to sample reference profiles in a grid
   */
  class GridSampler {
  public:

    /**
     * Initialise the sampler
     * @param image_size The size of the image to sample
     * @param ngrid The number of grid points
     */
    GridSampler(int2 image_size, int2 ngrid)
      : image_size_(image_size),
        ngrid_(ngrid),
        step_((double)image_size_[0] / (double)ngrid_[0],
              (double)image_size_[1] / (double)ngrid_[1]) {
      DIALS_ASSERT(image_size.all_gt(0));
      DIALS_ASSERT(ngrid.all_gt(0));
    }

    /**
     * @returns The image size of the grid.
     */
    int2 image_size() const {
      return image_size_;
    }

    /**
     * @returns the number of grid points.
     */
    int2 ngrid() const {
      return ngrid_;
    }

    /**
     * @returns The step size
     */
    double2 step() const {
      return step_;
    }

    /**
     * @returns The total number of grid points
     */
    std::size_t size() const {
      return ngrid_[0] * ngrid_[1];
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xy The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(double2 xy) const {
      int ix = (int)floor(xy[0] / step_[0]);
      int iy = (int)floor(xy[1] / step_[1]);
      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (ix >= image_size_[0]) ix = image_size_[0] - 1;
      if (iy >= image_size_[1]) iy = image_size_[1] - 1;
      return ix + iy * ngrid_[0];
    }

    /**
     * Get the x, y coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y coordinate of the profile
     */
    double2 operator[](std::size_t index) const {
      double x = ((index % ngrid_[0]) + 0.5) * step_[0];
      double y = ((index / ngrid_[0]) + 0.5) * step_[1];
      return double2(x, y);
    }

  private:
    int2 image_size_;
    int2 ngrid_;
    double2 step_;
  };

}} // namespace dials::algorithms


#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_H */
