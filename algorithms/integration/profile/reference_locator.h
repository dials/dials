/*
 * reference_locator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LOCATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LOCATOR_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::small;
  using scitbx::af::int3;
  using scitbx::af::double3;
  using scitbx::af::flex_double_const_ref;

  /**
   * A class to provide access to reference profiles by querying either
   * profile index or image volume coordinate.
   */
  template <typename ImageSampler>
  class ReferenceLocator {
  public:

    typedef ImageSampler sampler_type;

    /**
     * Instantiate the class with the reference prodiles and sampler
     * @param profiles The array of reference prodiles
     * @param sampler The sampler object.
     */
    ReferenceLocator(const flex_double_const_ref &profiles,
                     const ImageSampler &sampler)
      : profiles_(profiles),
        sampler_(sampler) {
      DIALS_ASSERT(profiles.accessor().all().size() == 4);
      DIALS_ASSERT(profiles.accessor().all().all_gt(0));
      DIALS_ASSERT(profiles.accessor().all()[0] == sampler_.size());
    }

    /** @returns The number of profiles. */
    std::size_t size() const {
      return sampler_.size();
    }

    /** @returns The sampler object */
    ImageSampler sampler() const {
      return sampler_;
    }

    /** @returns The whole list of profiles */
    flex_double_const_ref profile() const {
      return profiles_;
    }

    /**
     * Get the index of the profile nearest to the given coordinate.
     * @param xyz The image volume coordinate
     * @returns The index of the nearest profile to point xyz
     */
    std::size_t index(double3 xyz) const {
      return sampler_.nearest(xyz);
    }

    /**
     * Get the profile at the given index
     * @param index The profile index
     * @returns The profile array
     */
    flex_double_const_ref profile(std::size_t index) const {
      DIALS_ASSERT(index < sampler_.size());
      small<long, 10> profile_size = profiles_.accessor().all();
      DIALS_ASSERT(profile_size.size() == 4);
      std::size_t num = profile_size[1] * profile_size[2] * profile_size[3];
      return flex_double_const_ref(&profiles_[index], num);
    }

    /**
     * Get the profile at the given detector coordinate
     * @param xyz The detector coordinate
     * @returns The profile array
     */
    flex_double_const_ref profile(double3 xyz) const {
      std::size_t index = sampler_.nearest(xyz);
      return profile(index);
    }

    /**
     * Get the coordinate of the profile at the given index
     * @param index The profile index
     * @returns The profile coordinate
     */
    double3 coord(std::size_t index) const {
      return sampler_[index];
    }

    /**
     * Get the coordinate of the profile nearest the given coordinate
     * @param xyz The volume coordinate
     * @returns The coordinate of the profile.
     */
    double3 coord(double3 xyz) const {
      return sampler_[sampler_.nearest(xyz)];
    }

  private:
    flex_double_const_ref profiles_;
    ImageSampler sampler_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LOCATOR_H */
