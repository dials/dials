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
  template <typename FloatType, typename ImageSampler>
  class ReferenceLocator {
  public:

    typedef FloatType float_type;
    typedef ImageSampler sampler_type;

    /**
     * Instantiate the class with the reference prodiles and sampler
     * @param profiles The array of reference prodiles
     * @param sampler The sampler object.
     */
    ReferenceLocator(const af::versa< FloatType, af::c_grid<4> > &profiles,
                     const af::versa< bool, af::c_grid<4> > &masks,
                     const ImageSampler &sampler)
      : profiles_(profiles),
        masks_(masks),
        sampler_(sampler) {
      DIALS_ASSERT(profiles.accessor().all_eq(masks.accessor()));
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

    /**
     * Get the indices of the nearest profiles used in learning.
     * @param xyz The image volume coordinate
     * @return The list of indices
     */
    af::shared<std::size_t> indices(double3 xyz) const {
      return sampler_.nearest_n(xyz);
    }

    /** @returns The whole list of profiles */
    af::versa< FloatType, af::c_grid<4> > profile() const {
      return profiles_;
    }

    /** @returns The profile masks */
    af::versa< bool, af::c_grid<4> > mask() const {
      return masks_;
    }

    /**
     * Get the profile at the given index
     * @param index The profile index
     * @returns The profile array
     */
    af::const_ref< FloatType, af::c_grid<3> > profile_ref(std::size_t index) const {

      // Check the index and size of the profile array
      DIALS_ASSERT(index < sampler_.size());
      int4 profile_size = profiles_.accessor();

      // Return the result
      std::size_t j = index*profile_size[3]*profile_size[2]*profile_size[1];
      return af::const_ref< FloatType, af::c_grid<3> >(
          &profiles_[j],
          af::c_grid<3>(
            profile_size[1],
            profile_size[2],
            profile_size[3]));
    }

    af::const_ref< bool, af::c_grid<3> > mask_ref(std::size_t index) const {

      // Check the index and size of the profile array
      DIALS_ASSERT(index < sampler_.size());
      int4 mask_size = masks_.accessor();

      // Return the result
      std::size_t j = index*mask_size[3]*mask_size[2]*mask_size[1];
      return af::const_ref< bool, af::c_grid<3> >(
          &masks_[j],
          af::c_grid<3>(
            mask_size[1],
            mask_size[2],
            mask_size[3]));
    }

    /**
     * Get the profile at the given index
     * @param index The profile index
     * @returns The profile array
     */
    af::versa< FloatType, af::c_grid<3> > profile(std::size_t index) const {

      // Check the index and size of the profile array
      DIALS_ASSERT(index < sampler_.size());
      int4 profile_size = profiles_.accessor();

      // Unfortunately, you can't take a reference from a versa array and
      // return to python so we'll just have to make a copy.
      af::versa< FloatType, af::c_grid<3> > result(
        af::c_grid<3>(profile_size[1], profile_size[2], profile_size[3]),
        af::init_functor_null<FloatType>());
      std::size_t j = index*profile_size[3]*profile_size[2]*profile_size[1];
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = profiles_[j + i];
      }

      // Return the result
      return result;
    }

    /**
     * Get the mask at the given index
     * @param index The profile index
     * @returns The mask array
     */
    af::versa< bool, af::c_grid<3> > mask(std::size_t index) const {

      // Check the index and size of the profile array
      DIALS_ASSERT(index < sampler_.size());
      int4 mask_size = masks_.accessor();

      // Unfortunately, you can't take a reference from a versa array and
      // return to python so we'll just have to make a copy.
      af::versa< bool, af::c_grid<3> > result(
        af::c_grid<3>(mask_size[1], mask_size[2], mask_size[3]),
        af::init_functor_null<bool>());
      std::size_t j = index*mask_size[3]*mask_size[2]*mask_size[1];
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = masks_[j + i];
      }

      // Return the result
      return result;
    }

    /**
     * Get the profile at the given detector coordinate
     * @param xyz The detector coordinate
     * @returns The profile array
     */
    af::const_ref< FloatType, af::c_grid<3> > profile_ref(double3 xyz) const {
      return profile_ref(sampler_.nearest(xyz));
    }
    af::versa< FloatType, af::c_grid<3> > profile(double3 xyz) const {
      return profile(sampler_.nearest(xyz));
    }

    /**
     * Get the profile at the given detector coordinate
     * @param xyz The detector coordinate
     * @returns The profile array
     */
    af::const_ref< bool, af::c_grid<3> > mask_ref(double3 xyz) const {
      return mask_ref(sampler_.nearest(xyz));
    }
    af::versa< bool, af::c_grid<3> > mask(double3 xyz) const {
      return mask(sampler_.nearest(xyz));
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

    /**
     * Get the weight at a point for a reference profile
     * @param index The index of the profile
     * @patam coord The coordinate of the point
     */
    double weight(std::size_t index, double3 xyz) const {
      return sampler_.weight(index, xyz);
    }

    /**
     * @returns The correlations between the reference profiles
     */
    af::versa< double, af::c_grid<2> > correlations() const {
      af::versa< double, af::c_grid<2> > result(af::c_grid<2>(size(), size()));
      for (std::size_t j = 0; j < size(); ++j) {
        af::versa<FloatType, af::c_grid<3> > pj = profile(j);
        for (std::size_t i = 0; i < size(); ++i) {
          af::versa<FloatType, af::c_grid<3> > pi = profile(i);
          result(j,i) = compute_correlation(pj.const_ref(), pi.const_ref());
        }
      }
      return result;
    }

  private:

    /*
     * Compute the correlation coefficient between the profile and reference
     */
    double
    compute_correlation(const af::const_ref<FloatType, af::c_grid<3> > &x,
                        const af::const_ref<FloatType, af::c_grid<3> > &y) const {
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < x.size(); ++i) {
        xb += x[i];
        yb += y[i];
        count++;
      }
      DIALS_ASSERT(count > 0);
      xb /= count;
      yb /= count;
      double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
      for (std::size_t i = 0; i < x.size(); ++i) {
        double dx = x[i] - xb;
        double dy = y[i] - yb;
        sdxdy += dx*dy;
        sdx2 += dx*dx;
        sdy2 += dy*dy;
      }
      DIALS_ASSERT(sdx2 > 0.0 && sdy2 > 0.0);
      return sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
    }

    af::versa< FloatType, af::c_grid<4> > profiles_;
    af::versa< bool, af::c_grid<4> > masks_;
    ImageSampler sampler_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LOCATOR_H */
