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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using scitbx::af::double3;

  /**
   * Class to sample reference profiles in a grid
   */
  class GridSampler {
  public:

    /**
     * Initialise the sampler
     * @param volume_size The size of the volume to sample
     * @param grid_size The number of grid points
     */
    GridSampler(int3 volume_size, int3 grid_size)
      : volume_size_(volume_size),
        grid_size_(grid_size),
        step_size_((double)volume_size_[0] / (double)grid_size_[0],
                   (double)volume_size_[1] / (double)grid_size_[1],
                   (double)volume_size_[2] / (double)grid_size_[2]) {
      DIALS_ASSERT(volume_size.all_gt(0));
      DIALS_ASSERT(grid_size.all_gt(0));
    }

    /**
     * @returns The volume size of the grid.
     */
    int3 volume_size() const {
      return volume_size_;
    }

    /**
     * @returns the number of grid points.
     */
    int3 grid_size() const {
      return grid_size_;
    }

    /**
     * @returns The step size
     */
    double3 step_size() const {
      return step_size_;
    }

    /**
     * @returns The total number of grid points
     */
    std::size_t size() const {
      return grid_size_[0] * grid_size_[1] * grid_size_[2];
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xyz The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(double3 xyz) const {
      int ix = (int)floor(xyz[0] / step_size_[0]);
      int iy = (int)floor(xyz[1] / step_size_[1]);
      int iz = (int)floor(xyz[2] / step_size_[2]);
      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (iz < 0) iz = 0;
      if (ix >= grid_size_[0]) ix = grid_size_[0] - 1;
      if (iy >= grid_size_[1]) iy = grid_size_[1] - 1;
      if (iz >= grid_size_[2]) iz = grid_size_[2] - 1;
      return index(ix, iy, iz);
    }

    /**
     * Find the nearest n reference profiles to the given point.
     * This is mostly used during learning to get the neighbouring profiles.
     * @param xyz The coordinate
     * @returns A list of reference profile indices
     */
    af::shared<std::size_t> nearest_n(double3 xyz) const {
      int ix = (int)floor(xyz[0] / step_size_[0]);
      int iy = (int)floor(xyz[1] / step_size_[1]);
      int iz = (int)floor(xyz[2] / step_size_[2]);
      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (iz < 0) iz = 0;
      if (ix >= grid_size_[0]) ix = grid_size_[0] - 1;
      if (iy >= grid_size_[1]) iy = grid_size_[1] - 1;
      if (iz >= grid_size_[2]) iz = grid_size_[2] - 1;
      af::shared<std::size_t> result;
      result.push_back(index(ix, iy, iz));
      //if (ix > 0) result.push_back(index(ix-1, iy, iz));
      //if (iy > 0) result.push_back(index(ix, iy-1, iz));
      //if (iz > 0) result.push_back(index(ix, iy, iz-1));
      //if (ix < grid_size_[0]-1) result.push_back(index(ix+1, iy, iz));
      //if (iy < grid_size_[1]-1) result.push_back(index(ix, iy+1, iz));
      //if (iz < grid_size_[2]-1) result.push_back(index(ix, iy, iz+1));
      return result;
    }

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    double3 operator[](std::size_t index) const {
      DIALS_ASSERT(index < size());
      int i = index % grid_size_[0];
      int jk = index / grid_size_[0];
      int j = jk % grid_size_[1];
      int k = jk / grid_size_[1];
      double x = (i + 0.5) * step_size_[0];
      double y = (j + 0.5) * step_size_[1];
      double z = (k + 0.5) * step_size_[2];
      return double3(x, y, z);
    }

  private:

    /**
     * Create a profile index
     */
    std::size_t index(std::size_t ix, std::size_t iy, std::size_t iz) const {
      return ix + iy * grid_size_[0] + iz * grid_size_[0] * grid_size_[1];
    }

    int3 volume_size_;
    int3 grid_size_;
    double3 step_size_;
  };

}} // namespace dials::algorithms


#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_GRID_SAMPLER_H */
