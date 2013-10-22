/*
 * reference_learner.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H

#include <dials/model/data/reflection.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int3;
  using scitbx::af::int4;
  using dials::model::Reflection;

  /**
   * Class to learn the reference profiles
   */
  template <typename Sampler>
  class ReferenceLearner {
  public:

    typedef Reflection::float_type float_type;
    typedef Sampler sampler_type;
    typedef ReferenceLocator<float_type, sampler_type> locator_type;

    /**
     * Initialise the learner class
     * @param sampler The image volume sampler
     * @param grid_size The size of the grid
     * @param threshold The signal threshold
     */
    ReferenceLearner(const sampler_type &sampler,
                     int3 grid_size, double threshold)
      : locator_(allocate_profiles(sampler.size(), grid_size), sampler),
        threshold_(threshold) {}

    /**
     * @returns The reference profile locator.
     */
    locator_type locate() {
      return locator_;
    }

    /**
     * Learn the reference profiles from the reflection list.
     * @param reflections The list of reflections
     */
    void learn(const af::const_ref<Reflection> reflections) {
      // Add the contributions of all the reflections to the references
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].is_valid()) {
          add_reflection(reflections[i]);
        }
      }

      // Normalize the reference profiles
      normalize_reference_profiles();
    }

  private:

    /**
     * Allocate the array of reference profiles
     * @param num The number of profiles
     * @param grid_size The size of each profile
     * @returns The array of reference profiles
     */
    af::versa< float_type, af::c_grid<4> > allocate_profiles(
        std::size_t num, int3 grid_size) {
      DIALS_ASSERT(num > 0);
      DIALS_ASSERT(grid_size.all_gt(0));
      return af::versa< float_type, af::c_grid<4> >(af::c_grid<4>(
        int4(num, grid_size[0], grid_size[1], grid_size[2])), 0);
    }

    /**
     * Add a reflection to add to the reference profile learning.
     * @param reflection The reflection to add
     */
    void add_reflection(const Reflection &reflection) {
      vec2<double> image_coord = reflection.get_image_coord_px();
      double frame_number = reflection.get_frame_number();
      vec3<double> coord(image_coord[0], image_coord[1], frame_number);
      add_reflection(reflection.get_transformed_shoebox().const_ref(), coord);
    }

    /**
     * Add the actual reflection data to the profile
     * @param profile The reflection profile
     * @param coord The coordinate of the reflection
     */
    void add_reflection(const af::const_ref< float_type, af::c_grid<3> > profile,
        vec3<double> coord) {

      // Get the expected profile size
      int4 size_all = locator_.profile().accessor();
      int3 size(size_all[1], size_all[2], size_all[3]);

      // Ensure that the profiles are the correct size
      DIALS_ASSERT(profile.accessor().all_eq(size));

      // Find the nearest reference profile
      std::size_t index = locator_.index(coord);
      vec3<double> coord_b = locator_.coord(index);

      // Get the reference profile
      af::ref<float_type> reference = reference_profile(index);

      // Calculate the weighting by distance, ensure we don't get silly weights
      // for really close reflections by setting minimum distance to 1.
      double distance = (coord - coord_b).length();
      double weight = 1.0 / (distance < 1.0 ? 1.0 : distance);

      // Calculate the sum of intensity to normalize the reflection profile
      double sum_profile = 0.0;
      for (std::size_t i = 0; i < profile.size(); ++i) {
        sum_profile += profile[i];
      }

      // Add to the reference profile
      if (sum_profile > 0) {
        for (std::size_t i = 0; i < reference.size(); ++i) {
          reference[i] += weight * profile[i] / sum_profile;
        }
      }
    }

    /**
     * Normalize the reference profiles.
     */
    void normalize_reference_profiles() {
      // Normalize all the reference profiles
      for (std::size_t i = 0; i < locator_.size(); ++i) {
        normalize_reference_profile(i);
      }
    }

    /**
     * Normalize the reference profile at the given index
     * @param index The reference profile index
     */
    void normalize_reference_profile(std::size_t index) {

      // Get the reference profile at the index
      af::ref<float_type> reference = reference_profile(index);

      // Calculate the profile maximum and signal threshold
      double profile_maximum = max(reference);
      double threshold = threshold_ * profile_maximum;

      // Get the sum of signal pixels
      double signal_sum = 0.0;
      for (std::size_t i = 0; i < reference.size(); ++i) {
        if (reference[i] >= threshold) {
          signal_sum += reference[i];
        } else {
          reference[i] = 0.0;
        }
      }

      // If the reference profile sum is <= 0 then return
      DIALS_ASSERT(signal_sum > 0);

      // Normalize the profile such that sum of signal pixels == 1
      for (std::size_t i = 0; i < reference.size(); ++i) {
        reference[i] /= signal_sum;
      }
    }

    /**
     * Get a reference to the reflection profile
     * @param index The index of the profile to get
     * @returns The reference to the profile.
     */
    af::ref<float_type> reference_profile(std::size_t index) {
      af::versa<float_type, af::c_grid<4> > all_profiles = locator_.profile();
      int4 size = all_profiles.accessor();
      int offset = size[1] * size[2] * size[3];
      return af::ref<float_type>(&all_profiles[index * offset], offset);
    }

    locator_type locator_;
    double threshold_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H */
