/*
 * grid_sampler_2d.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_2D_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_2D_H

#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::double2;

  /**
   * Class to sample reference profiles in a grid
   */
  class GridSampler2D {
  public:

    /**
     * Initialise the sampler
     * @param image_size The size of the volume to sample
     * @param grid_size The number of grid points
     */
    GridSampler2D(int2 image_size, int2 grid_size)
      : image_size_(image_size),
        grid_size_(grid_size),
        step_size_((double)image_size_[0] / (double)grid_size_[0],
                   (double)image_size_[1] / (double)grid_size_[1]) {
      DIALS_ASSERT(image_size.all_gt(0));
      DIALS_ASSERT(grid_size.all_gt(0));
    }

    /**
     * @returns The volume size of the grid.
     */
    int2 image_size() const {
      return image_size_;
    }

    /**
     * @returns the number of grid points.
     */
    int2 grid_size() const {
      return grid_size_;
    }

    /**
     * @returns The step size
     */
    double2 step_size() const {
      return step_size_;
    }

    /**
     * @returns The total number of grid points
     */
    std::size_t size() const {
      return grid_size_[0] * grid_size_[1];
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xyz The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(double2 xy) const {
      DIALS_ASSERT(xy[0] >= 0 && xy[1] >= 0);
      DIALS_ASSERT(xy[0] < image_size_[0] && xy[1] < image_size_[1]);
      int ix = (int)floor(xy[0] / step_size_[0]);
      int iy = (int)floor(xy[1] / step_size_[1]);
      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (ix >= grid_size_[0]) ix = grid_size_[0] - 1;
      if (iy >= grid_size_[1]) iy = grid_size_[1] - 1;
      return index(ix, iy);
    }

    /**
     * Find the nearest n reference profiles to the given point.
     * This is mostly used during learning to get the neighbouring profiles.
     * @param xyz The coordinate
     * @returns A list of reference profile indices
     */
    af::shared<std::size_t> nearest_n(double2 xy) const {
      DIALS_ASSERT(xy[0] >= 0 && xy[1] >= 0);
      DIALS_ASSERT(xy[0] < image_size_[0] && xy[1] < image_size_[1]);
      double fx = xy[0] / step_size_[0];
      double fy = xy[1] / step_size_[1];
      int ix = (int)floor(fx);
      int iy = (int)floor(fy);
      DIALS_ASSERT(ix >= 0 && ix < grid_size_[0]);
      DIALS_ASSERT(iy >= 0 && iy < grid_size_[1]);
      double dx1 = std::abs(fx - ix);
      double dx2 = std::abs(fx - ix - 1);
      double dy1 = std::abs(fy - iy);
      double dy2 = std::abs(fy - iy - 1);
      int ix2 = dx1 < dx2 ? ix - 1 : ix + 1;
      int iy2 = dy1 < dy2 ? iy - 1 : iy + 1;
      bool xv = ix2 >= 0 && ix2 < grid_size_[0];
      bool yv = iy2 >= 0 && iy2 < grid_size_[1];
      af::shared<std::size_t> result;
      result.push_back(index(ix, iy));
      if (xv) {
        result.push_back(index(ix2, iy));
      }
      if (yv) {
        result.push_back(index(ix, iy2));
      }
      if (xv && yv) {
        result.push_back(index(ix2, iy2));
      }
      return result;
    }

    /**
     * Get the weight for the given profile at the given coordinate.
     * @param index The profile index
     * @param xyz The coordinate
     * @returns The weight (between 1.0 and 0.0)
     */
    double weight(std::size_t index, double2 xy) const {
      double2 c = (*this)[index];
      double dx = (c[0] - xy[0]) / step_size_[0];
      double dy = (c[1] - xy[1]) / step_size_[1];
      double d = std::sqrt(dx*dx + dy*dy);
      return std::exp(-4.0*d*d*std::log(2.0));
    }

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    double2 operator[](std::size_t index) const {
      DIALS_ASSERT(index < size());
      int i = index % grid_size_[0];
      int j = index / grid_size_[0];
      double x = (i + 0.5) * step_size_[0];
      double y = (j + 0.5) * step_size_[1];
      return double2(x, y);
    }

  private:

    /**
     * Create a profile index
     */
    std::size_t index(std::size_t ix, std::size_t iy) const {
      return ix + iy * grid_size_[0];
    }


    int2 image_size_;
    int2 grid_size_;
    double2 step_size_;
  };

}} // namespace dials::algorithms


#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_H */

