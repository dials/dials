/*
 * xds_circle_sampler.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_XDS_CIRCLE_SAMPLER_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_XDS_CIRCLE_SAMPLER_H

#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::constants::pi;
  using scitbx::constants::two_pi;
  using scitbx::af::int3;
  using scitbx::af::double2;
  using scitbx::af::double3;
  using scitbx::af::min;

  /**
   * Class to sample image with reference profiles as in XDS
   */
  class XdsCircleSampler {
  public:

    /**
     * Initialise the sampler
     * @param volume_size The size of the image to sample
     * @param num_z The number of grid points in z
     */
    XdsCircleSampler(int3 volume_size, std::size_t num_z)
      : volume_size_(volume_size),
        centre_(volume_size_[0] / 2.0, volume_size_[1] / 2.0),
        num_z_(num_z),
        nprofile_(9) {
      DIALS_ASSERT(volume_size_.all_gt(0));
      DIALS_ASSERT(num_z > 0);
      r0_ = min(centre_.const_ref());
      r1_ = r0_ / 3.0;
      r2_ = r1_ * std::sqrt(5.0);
      step_size_ = (double)volume_size_[2] / (double)num_z_;
    }

    /**
     * @returns The image size of the grid.
     */
    int3 volume_size() const {
      return volume_size_;
    }

    /**
     * @returns The number of grid points in z
     */
    std::size_t num_z() const {
      return num_z_;
    }

    /**
     * @returns The image centre.
     */
    double2 image_centre() const {
      return centre_;
    }

    /**
     * @returns The outer radius
     */
    double r0() const {
      return r0_;
    }

    /**
     * @returns The radius assigned to the inner profile
     */
    double r1() const {
      return r1_;
    }

    /**
     * @returns The radius of the outer reference profiles
     */
    double r2() const {
      return r2_;
    }

    /**
     * @returns The total number of reference profiles
     */
    std::size_t size() const {
      return nprofile_ * num_z_;
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xyz The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(double3 xyz) const {

      // Get the index at the image coordinate
      std::size_t ij = index_at_image_coord(double2(xyz[0], xyz[1]));

      // Get the k index
      int k = (int)floor(xyz[2] / step_size_);
      if (k < 0) k = 0;
      if (k >= num_z_) k = num_z_ - 1;

      // Return the index
      return ij + k * nprofile_;
    }

    /**
     * Get the x, y, z coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y, z coordinate of the profile
     */
    double3 operator[](std::size_t index) const {

      // Ensure index is within range
      DIALS_ASSERT(index >= 0 && index < size());

      // Get the image coordinate at the index
      double2 xy = image_coord_at_index(index % nprofile_);

      // Calculate the z coordinate
      double z = (index / nprofile_ + 0.5) * step_size_;

      // Return the x, y, z coordinate
      return double3(xy[0], xy[1], z);
    }

  private:

    /**
     * Find the nearest reference profile to the given point.
     * @param xy The coordinate
     * @returns The index of the reference profile
     */
    std::size_t index_at_image_coord(double2 xy) const {

      // Get the radius and angle at the point
      double xmc = xy[0] - centre_[0];
      double ymc = xy[1] - centre_[1];
      double r = std::sqrt(xmc*xmc + ymc*ymc);
      double t = atan2(ymc, xmc);

      // If radius is less than the inner radius return 0
      if (r < r1_) return 0;

      // Calculate the index of the profile about the circle
      int angular_index = (int)floor(t * (nprofile_ - 1) / two_pi + 0.5);

      // Return the index
      return (angular_index % (nprofile_ - 1)) + 1;
    }

    /**
     * Get the x, y coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y coordinate of the profile
     */
    double2 image_coord_at_index(std::size_t index) const {
      // Ensure index is within range
      DIALS_ASSERT(index >= 0 && index < nprofile_);

      // If index is zero then return central reference position
      if (index == 0) {
          return centre_;
      }

      // Calculate the angle and get the x, y coordinate of the profile
      double theta = (index - 1) * two_pi / (nprofile_ - 1);
      double x = centre_[0] + r2_ * cos(theta);
      double y = centre_[1] + r2_ * sin(theta);
      return double2(x, y);
    }

    int3 volume_size_;
    double2 centre_;
    std::size_t num_z_;
    std::size_t nprofile_;
    double step_size_;
    double r0_, r1_, r2_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_XDS_CIRCLE_SAMPLER_H */
