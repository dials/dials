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
  using scitbx::af::int2;
  using scitbx::af::double2;
  using scitbx::af::min;

  /**
   * Class to sample image with reference profiles as in XDS
   */
  class XdsCircleSampler {
  public:

    /**
     * Initialise the sampler
     * @param image_size The size of the image to sample
     */
    XdsCircleSampler(int2 image_size)
      : image_size_(image_size),
        centre_(image_size_[0] / 2.0, image_size_[1] / 2.0),
        nprofile_(9) {
      DIALS_ASSERT(image_size_.all_gt(0));
      r0_ = min(centre_.const_ref());
      r1_ = r0_ / 3.0;
      r2_ = r1_ * sqrt(5.0);
    }

    /**
     * @returns The image size of the grid.
     */
    int2 image_size() const {
      return image_size_;
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
      return nprofile_;
    }

    /**
     * Find the nearest reference profile to the given point.
     * @param xy The coordinate
     * @returns The index of the reference profile
     */
    std::size_t nearest(double2 xy) const {

      // Get the radius and angle at the point
      double xmc = xy[0] - centre_[0];
      double ymc = xy[1] - centre_[1];
      double r = sqrt(xmc*xmc + ymc*ymc);
      double t = atan2(ymc, xmc) + pi;

      // If radius is less than the inner radius return 0
      if (r < r1_) return 0;

      // Calculate the index of the profile about the circle
      int angular_index = (int)floor(t * (nprofile_ - 1) / two_pi + 0.5);

      // Return the index
      return angular_index % (nprofile_ - 1) + 1;

    }

    /**
     * Get the x, y coordinate of the reference profile at the given index.
     * @param index The index of the reference profile.
     * @returns The x, y coordinate of the profile
     */
    double2 operator[](std::size_t index) const {

      // Ensure index is within range
      DIALS_ASSERT(index >= 0 && index < nprofile_);

      // If index is zero then return central reference position
      if (index == 0) {
          return centre_;
      }

      // Calculate the angle and get the x, y coordinate of the profile
      double theta = (index - 1) * two_pi / (nprofile_ - 1);
      double x = centre_[0] + r2_ * sin(theta);
      double y = centre_[1] + r2_ * cos(theta);
      return double2(x, y);
    }

  private:
    int2 image_size_;
    double2 centre_;
    std::size_t nprofile_;
    double r0_, r1_, r2_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_XDS_CIRCLE_SAMPLER_H */
