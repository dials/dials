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

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <dials/model/data/transformed_shoebox.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int3;
  using scitbx::af::int4;
  using model::TransformedShoebox;

  /**
   * Class to learn the reference profiles
   */
  template <typename Sampler>
  class ReferenceLearner {
  public:

    typedef double float_type;
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
      : locator_(
          allocate_profiles<float_type>(sampler.size(), grid_size, 0.0),
          allocate_profiles<bool>(sampler.size(), grid_size, true),
          sampler),
        threshold_(threshold),
        counts_(sampler.size(), 0) {}

    /**
     * @returns The reference profile locator.
     */
    locator_type locate() {
      return locator_;
    }

    /**
     * Learn the reference profiles from the reflection list.
     * @param profiles The list of profiles
     */
    void learn(const af::const_ref< TransformedShoebox > profiles,
               const af::const_ref< vec3<double> > coords) {
      // Add the contributions of all the profiles to the references
      DIALS_ASSERT(profiles.size() == coords.size());
      for (std::size_t i = 0; i < profiles.size(); ++i) {
        add_profile(profiles[i], coords[i]);
      }

      // Normalize the reference profiles
      normalize_reference_profiles();
    }

    /**
     * @returns The number of reflections contributing to each profile.
     */
    af::shared<int> counts() const {
      af::shared<int> result(counts_.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = counts_[i];
      }
      return result;
    }

  private:

    /**
     * Allocate the array of reference profiles
     * @param num The number of profiles
     * @param grid_size The size of each profile
     * @returns The array of reference profiles
     */
    template <typename T>
    af::versa< T, af::c_grid<4> > allocate_profiles(
        std::size_t num, int3 grid_size, T value) {
      DIALS_ASSERT(num > 0);
      DIALS_ASSERT(grid_size.all_gt(0));
      return af::versa< T, af::c_grid<4> >(af::c_grid<4>(
        int4(num, grid_size[0], grid_size[1], grid_size[2])), value);
    }

    /**
     * Add a reflection to add to the reference profile learning.
     * @param profile The profile to add
     */
    void add_profile(const TransformedShoebox &profile, vec3<double> coord) {
      add_profile(profile.data.const_ref(), coord);
    }

    /**
     * Add the actual reflection data to the profile
     * @param profile The reflection profile
     * @param coord The coordinate of the reflection
     */
    void add_profile(const af::const_ref< float_type, af::c_grid<3> > profile,
        vec3<double> coord) {

      // Get the expected profile size
      int4 size_all = locator_.profile().accessor();
      int3 size(size_all[1], size_all[2], size_all[3]);

      // Ensure that the profiles are the correct size
      DIALS_ASSERT(profile.accessor().all_eq(size));

      // Calculate the sum of intensity to normalize the reflection profile
      double sum_profile = 0.0;
      for (std::size_t i = 0; i < profile.size(); ++i) {
        sum_profile += profile[i];
      }

      // Find the nearest reference profile
      if (sum_profile > 0) {
        af::shared<std::size_t> indices = locator_.indices(coord);
        for (std::size_t ii = 0; ii < indices.size(); ++ii) {
          std::size_t index = indices[ii];
          //vec3<double> coord_b = locator_.coord(index);
          counts_[index]++;

          // Get the reference profile
          af::ref<float_type> reference = reference_profile(index);

          // Calculate the weighting by distance, ensure we don't get silly weights
          // for really close reflections by setting minimum distance to 1.
          //double distance = (coord - coord_b).length();
          double weight = locator_.weight(index, coord);
          //double weight = (distance < 1.0 ? 1.0 : 1.0 / distance);

          // Add to the reference profile
          for (std::size_t i = 0; i < reference.size(); ++i) {
            reference[i] += weight * profile[i] / sum_profile;
          }
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
      af::ref<bool> mask = reference_mask(index);

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
          mask[i] = false;
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

    /**
     * Get a reference to the reflection profile
     * @param index The index of the profile to get
     * @returns The reference to the profile.
     */

    af::ref<bool> reference_mask(std::size_t index) {
      af::versa<bool, af::c_grid<4> > all_masks = locator_.mask();
      int4 size = all_masks.accessor();
      int offset = size[1] * size[2] * size[3];
      return af::ref<bool>(&all_masks[index * offset], offset);
    }

    locator_type locator_;
    double threshold_;
    af::shared<int> counts_;
  };

  /**
   * Class to learn the reference profiles
   */
  template <typename Sampler>
  class ReferenceLearnerNew {
  public:

    typedef double float_type;
    typedef Sampler sampler_type;
    typedef ReferenceLocator<float_type, sampler_type> locator_type;

    /**
     * Initialise the learner class
     * @param sampler The image volume sampler
     * @param grid_size The size of the grid
     * @param threshold The signal threshold
     */
    ReferenceLearnerNew(const sampler_type &sampler,
                     int3 grid_size, double threshold)
      : locator_(
          allocate_profiles<float_type>(sampler.size(), grid_size, 0.0),
          allocate_profiles<bool>(sampler.size(), grid_size, true),
          sampler),
        threshold_(threshold),
        counts_(sampler.size(), 0),
        finalized_(false) {}

    /**
     * @returns The reference profile locator.
     */
    locator_type locate() const {
      DIALS_ASSERT(finalized_);
      return locator_;
    }

    /**
     * Learn the reference profiles from the reflection list.
     * @param profiles The list of profiles
     */
    void add(const af::const_ref< float_type, af::c_grid<3> > profile,
        vec3<double> coord) {
      add_profile(profile, coord);
    }

    void finalize() {
      normalize_reference_profiles();
      finalized_ = true;
    }

    /**
     * @returns The number of reflections contributing to each profile.
     */
    af::shared<int> counts() const {
      af::shared<int> result(counts_.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = counts_[i];
      }
      return result;
    }

  private:

    /**
     * Allocate the array of reference profiles
     * @param num The number of profiles
     * @param grid_size The size of each profile
     * @returns The array of reference profiles
     */
    template <typename T>
    af::versa< T, af::c_grid<4> > allocate_profiles(
        std::size_t num, int3 grid_size, T value) {
      DIALS_ASSERT(num > 0);
      DIALS_ASSERT(grid_size.all_gt(0));
      return af::versa< T, af::c_grid<4> >(af::c_grid<4>(
        int4(num, grid_size[0], grid_size[1], grid_size[2])), value);
    }

    /**
     * Add the actual reflection data to the profile
     * @param profile The reflection profile
     * @param coord The coordinate of the reflection
     */
    void add_profile(const af::const_ref< float_type, af::c_grid<3> > profile,
        vec3<double> coord) {

      // Get the expected profile size
      int4 size_all = locator_.profile().accessor();
      int3 size(size_all[1], size_all[2], size_all[3]);

      // Ensure that the profiles are the correct size
      DIALS_ASSERT(profile.accessor().all_eq(size));

      // Calculate the sum of intensity to normalize the reflection profile
      double sum_profile = 0.0;
      for (std::size_t i = 0; i < profile.size(); ++i) {
        sum_profile += profile[i];
      }

      // Find the nearest reference profile
      if (sum_profile > 0) {
        af::shared<std::size_t> indices = locator_.indices(coord);
        for (std::size_t ii = 0; ii < indices.size(); ++ii) {
          std::size_t index = indices[ii];
          //vec3<double> coord_b = locator_.coord(index);
          counts_[index]++;

          // Get the reference profile
          af::ref<float_type> reference = reference_profile(index);

          // Calculate the weighting by distance, ensure we don't get silly weights
          // for really close reflections by setting minimum distance to 1.
          //double distance = (coord - coord_b).length();
          double weight = locator_.weight(index, coord);
          //double weight = (distance < 1.0 ? 1.0 : 1.0 / distance);

          // Add to the reference profile
          for (std::size_t i = 0; i < reference.size(); ++i) {
            reference[i] += weight * profile[i] / sum_profile;
          }
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
      af::ref<bool> mask = reference_mask(index);

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
          mask[i] = false;
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

    /**
     * Get a reference to the reflection profile
     * @param index The index of the profile to get
     * @returns The reference to the profile.
     */

    af::ref<bool> reference_mask(std::size_t index) {
      af::versa<bool, af::c_grid<4> > all_masks = locator_.mask();
      int4 size = all_masks.accessor();
      int offset = size[1] * size[2] * size[3];
      return af::ref<bool>(&all_masks[index * offset], offset);
    }

    locator_type locator_;
    double threshold_;
    af::shared<int> counts_;
    bool finalized_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H */
