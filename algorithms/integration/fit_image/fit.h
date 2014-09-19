/*
 * fit.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H
#define DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H

#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/integration/profile/fitting.h>


namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::af::ReferenceSpot;
  using dials::af::DontIntegrate;
  using dials::af::IntegratedPrf;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;

  /**
   * A class to specify the experiments input
   */
  class Spec {
  public:

    Spec(const Beam &beam,
         const Detector &detector,
         const Goniometer &goniometer,
         const Scan &scan,
         double delta_b,
         double delta_m)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        delta_b_(delta_b),
        delta_m_(delta_m) {
      DIALS_ASSERT(delta_b_ > 0);
      DIALS_ASSERT(delta_m_ > 0);
    }

    const Beam& beam() const {
      return beam_;
    }

    const Detector& detector() const {
      return detector_;
    }

    const Goniometer& goniometer() const {
      return goniometer_;
    }

    const Scan& scan() const {
      return scan_;
    }

    const double delta_b() const {
      return delta_b_;
    }

    const double delta_m() const {
      return delta_m_;
    }

  private:

    Beam beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    double delta_b_;
    double delta_m_;
  };


  /**
   * A class to compute a single reference profile.
   */
  template <typename FloatType>
  class SingleProfileLearner {
  public:

    typedef FloatType float_type;
    typedef af::versa< FloatType, af::c_grid<3> > profile_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;

    SingleProfileLearner(const Spec spec,
                         int3 grid_size)
      : spec_(spec),
        reference_data_(af::c_grid<3>(
              2 * grid_size[0] + 1,
              2 * grid_size[1] + 1,
              2 * grid_size[2] + 1)),
        reference_mask_(af::c_grid<3>(
              2 * grid_size[0] + 1,
              2 * grid_size[1] + 1,
              2 * grid_size[2] + 1)) {}

    void add(const Shoebox<> &sbox, const vec3<double> &xyz) {

    }

    profile_type get(
        std::size_t panel,
        vec3<double> xyz,
        int6 bbox) const {
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
      std::size_t xsize = bbox[1] - bbox[0];
      std::size_t ysize = bbox[3] - bbox[2];
      std::size_t zsize = bbox[5] - bbox[4];
      af::c_grid<3> accessor(zsize, ysize, xsize);
      return profile_type(accessor, 1);
    }

  private:

    Spec spec_;
    profile_type reference_data_;
    mask_type reference_mask_;
  };


  /**
   * A class to compute a single reference profile per experiment.
   */
  template <typename FloatType>
  class MultiExpSingleProfileLearner {
  public:

    typedef FloatType float_type;
    typedef SingleProfileLearner<FloatType> learner_type;
    typedef typename learner_type::profile_type profile_type;

    /**
     * Compute the reference profiles.
     * @param spec The list of experiment specifications
     * @param data The list of reflections.
     */
    MultiExpSingleProfileLearner(
        af::const_ref<Spec> spec,
        af::reflection_table data,
        std::size_t grid_size) {

      // Check the input
      DIALS_ASSERT(spec.size() > 0);
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(grid_size > 0);

      // Create the array of learners
      learner_.reserve(spec.size());
      for (std::size_t i = 0; i < spec.size(); ++i) {
        learner_.push_back(learner_type(spec[i],
              int3(grid_size, grid_size, grid_size)));
      }

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > xyzcal  = data["xyzcal.px"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Loop through all the reflections
      for (std::size_t i = 0; i < data.size(); ++i) {
        DIALS_ASSERT(id[i] < learner_.size());
        if (flags[i] & ReferenceSpot) {
          learner_[id[i]].add(shoebox[i], xyzcal[i]);
        }
      }
    }

    /**
     * Get the reference profile for the reflection.
     * @param id The experiment ID
     * @param panel The panel number
     * @param xyz The calculated xyz pixel coordinate
     * @param bbox The bounding box
     * @returns The reference profile for the reflection.
     */
    profile_type get(
        std::size_t id,
        std::size_t panel,
        vec3<double> xyz,
        int6 bbox) const {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id].get(panel, xyz, bbox);
    }

  private:

    std::vector<learner_type> learner_;
  };


  /**
   * A class to perform image space profile fitting.
   */
  class ImageSpaceProfileFitting {
  public:

    typedef Shoebox<>::float_type float_type;
    typedef MultiExpSingleProfileLearner<float_type> reference_learner_type;
    typedef reference_learner_type::profile_type profile_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;
    typedef ProfileFitting<float_type> fitting_type;

    ImageSpaceProfileFitting(std::size_t grid_size)
        : grid_size_(grid_size) {
      DIALS_ASSERT(grid_size > 0);
    }

    /**
     * Add some experiment data
     * @param spec The experiment specification
     */
    void add(const Spec &spec) {
      spec_.push_back(spec);
    }

    /**
     * Perform the profile fitting on the input data
     * @param data The reflection data
     */
    void execute(af::reflection_table data) const {

      // Check we have experimental setup data
      DIALS_ASSERT(spec_.size() > 0);

      // Check the input contains expected fields
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("xyzcal.px"));
      DIALS_ASSERT(data.contains("flags"));

      // Compute the reference profile
      reference_learner_type reference(spec_.const_ref(), data, grid_size_);

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > xyzcal  = data["xyzcal.px"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Get the new columns to set
      af::ref<double> intensity   = data["intensity.prf.value"];
      af::ref<double> variance    = data["intensity.prf.variance"];
      af::ref<double> correlation = data["profile_correlation"];

      // Do the profile fitting for all reflections
      for (std::size_t i = 0; i < data.size(); ++i) {

        // Initialise some stuff
        intensity[i] = 0;
        variance[i] = 0;
        correlation[i] = 0;
        flags[i] &= ~IntegratedPrf;

        // Integrate
        if (!(flags[i] & DontIntegrate)) {

          // Get the shoebox
          const Shoebox<> &sbox = shoebox[i];
          DIALS_ASSERT(sbox.is_consistent());

          // Get the profile for a given reflection
          profile_type profile = reference.get(
              id[i],
              sbox.panel,
              xyzcal[i],
              sbox.bbox);

          // Compute the integration mask
          mask_type mask(sbox.mask.accessor(), false);
          std::size_t mask_code = Valid | Foreground;
          std::size_t mask_count = 0;
          DIALS_ASSERT(profile.accessor().all_eq(sbox.mask.accessor()));
          for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
            if ((sbox.mask[j] & mask_code) == mask_code && profile[j] >= 0.0) {
              mask[j] = true;
              mask_count++;
            }
          }

          // Perform the profile fit
          if (mask_count > 0) {
            fitting_type fit(
                profile.const_ref(),
                mask.const_ref(),
                sbox.data.const_ref(),
                sbox.background.const_ref());

            // Set the data in the reflection
            intensity[i]   = fit.intensity();
            variance[i]    = fit.variance();
            correlation[i] = fit.correlation();

            // Set the integrated flag
            flags[i] |= IntegratedPrf;
          }
        }
      }
    }

  private:

    af::shared<Spec> spec_;
    std::size_t grid_size_;
  };

}} // namespace dials::algorithms


#endif // DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H
