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

#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using scitbx::af::int4;
  using scitbx::af::double3;

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
    ReferenceLocator(const af::versa< double, af::c_grid<4> > &profiles,
                     const ImageSampler &sampler)
      : profiles_(profiles),
        sampler_(sampler) {
      DIALS_ASSERT(profiles.accessor().all_gt(0));
      DIALS_ASSERT(profiles.accessor()[0] == sampler_.size());
    }

    /** @returns The number of profiles. */
    std::size_t size() const {
      return sampler_.size();
    }

    /** @returns The sampler object */
    ImageSampler sampler() const {
      return sampler_;
    }

    /**
     * Get the index of the profile nearest to the given coordinate.
     * @param xyz The image volume coordinate
     * @returns The index of the nearest profile to point xyz
     */
    std::size_t index(double3 xyz) const {
      return sampler_.nearest(xyz);
    }

    /** @returns The whole list of profiles */
    af::versa< double, af::c_grid<4> > profile() const {
      return profiles_;
    }

    /**
     * Get the profile at the given index
     * @param index The profile index
     * @returns The profile array
     */
    af::versa< double, af::c_grid<3> > profile(std::size_t index) const {

      // Check the index and size of the profile array
      DIALS_ASSERT(index < sampler_.size());
      int4 profile_size = profiles_.accessor();

      // Unfortunately, you can't take a reference from a versa array and
      // return to python so we'll just have to make a copy.
      af::versa< double, af::c_grid<3> > result(
        af::c_grid<3>(profile_size[1], profile_size[2], profile_size[2]),
        af::init_functor_null<double>());
      std::size_t j = index*profile_size[3]*profile_size[2]*profile_size[1];
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = profiles_[j + i];
      }

      // Return the result
      return result;
    }

    /**
     * Get the profile at the given detector coordinate
     * @param xyz The detector coordinate
     * @returns The profile array
     */
    af::versa< double, af::c_grid<3> > profile(double3 xyz) const {
      return profile(sampler_.nearest(xyz));
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
    af::versa< double, af::c_grid<4> > profiles_;
    ImageSampler sampler_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LOCATOR_H */
