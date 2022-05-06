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
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_GRID_SAMPLER_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_GRID_SAMPLER_H

#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <dials/algorithms/profile_model/modeller/sampler_interface.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::double3;
  using scitbx::af::int2;
  using scitbx::af::int3;

  /**
   * Class to sample reference profiles in a grid
   */
  class GridSampler : public SamplerIface {
  public:
    /**
     * Initialise the sampler
     * @param volume_size The size of the volume to sample
     * @param grid_size The number of grid points
     */
    GridSampler(int2 image_size, int2 scan_range, int3 grid_size)
        : image_size_(image_size),
          scan_range_(scan_range),
          scan_size_(scan_range[1] - scan_range[0]),
          grid_size_(grid_size) {
      // Check some input
      DIALS_ASSERT(image_size_.all_gt(0));
      DIALS_ASSERT(grid_size_.all_gt(0));
      DIALS_ASSERT(scan_range_[1] > scan_range_[0]);
      DIALS_ASSERT(scan_size_ > 0);

      // Compute the step size
      step_size_[0] = (double)image_size_[0] / (double)grid_size_[0];
      step_size_[1] = (double)image_size_[1] / (double)grid_size_[1];
      step_size_[2] = (double)scan_size_ / (double)grid_size_[2];
    }

    /**
     * @returns The volume size of the grid.
     */
    int2 image_size() const {
      return image_size_;
    }

    /**
     * @returns The scan range
     */
    int2 scan_range() const {
      return scan_range_;
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
    std::size_t nearest(std::size_t panel, double3 xyz) const {
      DIALS_ASSERT(panel == 0);
      DIALS_ASSERT(xyz[0] >= 0 && xyz[1] >= 0);
      DIALS_ASSERT(xyz[0] < image_size_[0]);
      DIALS_ASSERT(xyz[1] < image_size_[1]);
      DIALS_ASSERT(xyz[2] >= scan_range_[0]);
      DIALS_ASSERT(xyz[2] < scan_range_[1]);
      xyz[2] -= scan_range_[0];
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
    af::shared<std::size_t> nearest_n(std::size_t panel, double3 xyz) const {
      DIALS_ASSERT(panel == 0);
      std::size_t index = nearest(panel, xyz);
      af::shared<std::size_t> result = neighbours(index);
      result.push_back(index);
      return result;
    }

    /**
     * Get the weight for the given profile at the given coordinate.
     * @param index The profile index
     * @param xyz The coordinate
     * @returns The weight (between 1.0 and 0.0)
     */
    double weight(std::size_t index, std::size_t panel, double3 xyz) const {
      DIALS_ASSERT(panel == 0);
      double3 c = coord(index);
      double dx = (c[0] - xyz[0]) / step_size_[0];
      double dy = (c[1] - xyz[1]) / step_size_[1];
      double dz = (c[2] - xyz[2]) / step_size_[2];
      double d = std::sqrt(dx * dx + dy * dy + dz * dz);
      return std::exp(-4.0 * d * d * std::log(2.0));
    }

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    double3 coord(std::size_t index) const {
      DIALS_ASSERT(index < size());
      int i = index % grid_size_[0];
      int jk = index / grid_size_[0];
      int j = jk % grid_size_[1];
      int k = jk / grid_size_[1];
      double x = (i + 0.5) * step_size_[0];
      double y = (j + 0.5) * step_size_[1];
      double z = (k + 0.5) * step_size_[2] + scan_range_[0];
      return double3(x, y, z);
    }

    /**
     * Return the neighbouring grid points.
     */
    af::shared<std::size_t> neighbours(std::size_t index) const {
      DIALS_ASSERT(index < size());
      int i = index % grid_size_[0];
      int jk = index / grid_size_[0];
      int j = jk % grid_size_[1];
      int k = jk / grid_size_[1];
      int xs = grid_size_[0];
      int ys = grid_size_[1];
      int zs = grid_size_[2];
      af::shared<std::size_t> result;
      for (int kk = -1; kk <= 1; ++kk) {
        for (int jj = -1; jj <= 1; ++jj) {
          for (int ii = -1; ii <= 1; ++ii) {
            int iii = i + ii;
            int jjj = j + jj;
            int kkk = k + kk;
            if (iii >= 0 && iii < xs && jjj >= 0 && jjj < ys && kkk >= 0 && kkk < zs
                && !(kk == 0 && jj == 0 && ii == 0)) {
              result.push_back(iii + jjj * xs + kkk * xs * ys);
            }
          }
        }
      }
      return result;
    }

  private:
    /**
     * Create a profile index
     */
    std::size_t index(std::size_t ix, std::size_t iy, std::size_t iz) const {
      return ix + iy * grid_size_[0] + iz * grid_size_[0] * grid_size_[1];
    }

    int2 image_size_;
    int2 scan_range_;
    int scan_size_;
    int3 grid_size_;
    double3 step_size_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_GRID_SAMPLER_H */
