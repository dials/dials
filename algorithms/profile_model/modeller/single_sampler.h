/*
 * single_sampler.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SINGLE_SAMPLER_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SINGLE_SAMPLER_H

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
   * Class to sample a single reference profile per block
   */
  class SingleSampler : public SamplerIface {
  public:
    /**
     * Initialise the sampler
     * @param scan_range The scan range
     */
    SingleSampler(int2 scan_range, std::size_t grid_size)
        : scan_range_(scan_range),
          scan_size_(scan_range[1] - scan_range[0]),
          grid_size_(grid_size) {
      // Check some input
      DIALS_ASSERT(scan_range_[1] > scan_range_[0]);
      DIALS_ASSERT(scan_size_ > 0);
      DIALS_ASSERT(grid_size_ > 0);

      // Compute the step size
      step_size_ = (double)scan_size_ / (double)grid_size_;
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
    std::size_t grid_size() const {
      return grid_size_;
    }

    /**
     * @returns The step size
     */
    double step_size() const {
      return step_size_;
    }

    /**
     * @returns The total number of grid points
     */
    std::size_t size() const {
      return grid_size_;
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xyz The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(std::size_t panel, double3 xyz) const {
      DIALS_ASSERT(xyz[2] >= scan_range_[0]);
      DIALS_ASSERT(xyz[2] < scan_range_[1]);
      xyz[2] -= scan_range_[0];
      int iz = (int)floor(xyz[2] / step_size_);
      if (iz < 0) iz = 0;
      if (iz >= grid_size_) iz = grid_size_ - 1;
      return iz;
    }

    /**
     * Find the nearest n reference profiles to the given point.
     * This is mostly used during learning to get the neighbouring profiles.
     * @param xyz The coordinate
     * @returns A list of reference profile indices
     */
    af::shared<std::size_t> nearest_n(std::size_t panel, double3 xyz) const {
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
      double3 c = coord(index);
      double d = std::abs((c[2] - xyz[2]) / step_size_);
      return std::exp(-4.0 * d * d * std::log(2.0));
    }

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    double3 coord(std::size_t index) const {
      DIALS_ASSERT(index < size());
      double x = 0;
      double y = 0;
      double z = (index + 0.5) * step_size_ + scan_range_[0];
      return double3(x, y, z);
    }

    /**
     * Return the neighbouring grid points.
     */
    af::shared<std::size_t> neighbours(std::size_t index) const {
      DIALS_ASSERT(index < size());
      af::shared<std::size_t> result;
      if (index > 0) {
        result.push_back(index - 1);
      }
      if (index < grid_size_ - 1) {
        result.push_back(index + 1);
      }
      return result;
    }

  private:
    int2 scan_range_;
    int scan_size_;
    std::size_t grid_size_;
    double step_size_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_SINGLE_SAMPLER_H */
