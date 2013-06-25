#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H

#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int3;
  using scitbx::af::flex_double;
  using scitbx::af::flex_double_ref;
  using dials::model::ReflectionList;
  using dials::model::Reflection;

  template <typename Sampler>
  class ReferenceLearner {
  public:

    typedef Sampler sampler_type;
    typedef ReferenceLocator<sampler_type> locator_type;

    ReferenceLearner(const sampler_type &sampler, int3 grid_size, double threshold)
      : locator_(allocate_profiles(sampler.size(), grid_size), sampler),
        threshold_(threshold) {}

    locator_type locate() {
      return locator_;
    }

    void learn(const ReflectionList &reflections) {
      // Add the contributions of all the reflections to the references
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].get_status() == 0) {
          add_reflection(reflections[i]);
        }
      }

      // Normalize the reference profiles
      normalize_reference_profiles();
    }

  private:

    flex_double allocate_profiles(std::size_t num, int3 grid_size) {
      return flex_double(flex_grid<>(num, grid_size[0],
        grid_size[1], grid_size[2]), 0);
    }

    /**
     * Add a reflection to add to the reference profile learning.
     * @param reflection The reflection to add
     */
    void add_reflection(const Reflection &reflection) {
      vec2<double> image_coord = reflection.get_image_coord_px();
      double frame_number = reflection.get_frame_number();
      vec3<double> coord(image_coord[0], image_coord[1], frame_number);
      add_reflection(reflection.get_transformed_shoebox(), coord);
    }

    void add_reflection(flex_double profile, vec3<double> coord) {

      // Get the expected profile size
      small<long,10> size_all = locator_.profile().accessor().all();
      DIALS_ASSERT(size_all.size() == 4);
      small<long,10> size;
      size[0] = size_all[1];
      size[1] = size_all[2];
      size[2] = size_all[3];

      // Ensure that the profiles are the correct size
      DIALS_ASSERT(profile.accessor().all().size() == 3);
      DIALS_ASSERT(profile.accessor().all().all_eq(size));

      // Find the nearest reference profile
      std::size_t index = locator_.index(coord);
      vec3<double> coord_b = locator_.coord(index);

      // Get the reference profile
      flex_double_ref reference = reference_profile(index);

      // Calculate the weighting by distance
      double weight = 1.0 / (coord - coord_b).length();

      // Calculate the sum of intensity to normalize the reflection profile
      double sum_profile = 0.0;
      for (std::size_t i = 0; i < profile.size(); ++i) {
        sum_profile += profile[i];
      }

      // Add to the reference profile
      for (std::size_t i = 0; i < reference.size(); ++i) {
        reference[i] += weight * profile[i] / sum_profile;
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
      flex_double_ref reference = reference_profile(index);

      // Calculate the profile maximum and signal threshold
      double profile_maximum = max(reference);
      double threshold = threshold_ * profile_maximum;

      // Get the sum of signal pixels
      double signal_sum = 0.0;
      for (std::size_t i = 0; i < reference.size(); ++i) {
        if (reference[i] >= threshold) {
          signal_sum += reference[i];
        }
      }

      // Normalize the profile such that sum of signal pixels == 1
      for (std::size_t i = 0; i < reference.size(); ++i) {
        reference[i] /= signal_sum;
      }
    }

    flex_double_ref reference_profile(std::size_t index) {
      flex_double all_profiles = locator_.profile();
      small<long,10> size = all_profiles.accessor().all();
      DIALS_ASSERT(size.size() == 4);
      int offset = size[1] * size[2] * size[3];
      return flex_double_ref(&all_profiles[index * offset], offset);
    }

    locator_type locator_;
    double threshold_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_H */
