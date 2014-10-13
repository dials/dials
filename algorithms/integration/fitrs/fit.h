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

#ifndef DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H
#define DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H

#include <iostream>
#include <iomanip>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/integration/profile/fitting.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/polygon/clip/clip.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/algorithms/polygon/area.h>


namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dxtbx::model::Panel;
  using dials::af::ReferenceSpot;
  using dials::af::DontIntegrate;
  using dials::af::IntegratedPrf;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::algorithms::reflection_basis::CoordinateSystem;
  using dials::algorithms::polygon::simple_area;
  using dials::algorithms::polygon::clip::vert4;
  using dials::algorithms::polygon::clip::vert8;
  using dials::algorithms::polygon::clip::quad_with_convex_quad;
  using dials::algorithms::polygon::spatial_interpolation::reverse_quad_inplace_if_backward;

  template <typename T>
  T min4(T a, T b, T c, T d) {
    return std::min(std::min(a,b), std::min(c,d));
  }

  template <typename T>
  T max4(T a, T b, T c, T d) {
    return std::max(std::max(a,b), std::max(c,d));
  }

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
   * A class to perform reciprocal space profile fitting.
   */
  class ReciprocalSpaceProfileFitting {
  public:

    typedef Shoebox<>::float_type float_type;
    /* typedef MultiExpSingleProfileLearner<float_type> reference_learner_type; */
    /* typedef reference_learner_type::profile_type profile_type; */
    typedef af::versa< bool, af::c_grid<3> > mask_type;
    typedef ProfileFitting<float_type> fitting_type;

    ReciprocalSpaceProfileFitting(std::size_t grid_size)
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
    void
    execute(af::reflection_table data) const {

      // Check we have experimental setup data
      DIALS_ASSERT(spec_.size() > 0);

      // Check the input contains expected fields
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("s1"));
      DIALS_ASSERT(data.contains("xyzcal.mm"));
      DIALS_ASSERT(data.contains("flags"));

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > s1      = data["s1"];
      af::const_ref< vec3<double> > xyzcal  = data["xyzcal.mm"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Remove any reflections which have invalid pixels or are not fully
      // recorded from integration.
      for (std::size_t i = 0; i < data.size(); ++i) {
        Shoebox<> sbox = shoebox[i];
        std::size_t panel = sbox.panel;
        int6 bbox = sbox.bbox;
        DIALS_ASSERT(id[i] < spec_.size());
        int2 image_size = spec_[id[i]].detector()[panel].get_image_size();
        DIALS_ASSERT(sbox.is_consistent());
        if (bbox[0] < 0 || bbox[1] > image_size[0] ||
            bbox[2] < 0 || bbox[3] > image_size[1]) {
          flags[i] |= DontIntegrate;
          flags[i] &= ~ReferenceSpot;
          continue;
        }
        for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
          if (sbox.mask[j] & Foreground && !(sbox.mask[j] & Valid)) {
            flags[i] |= DontIntegrate;
            flags[i] &= ~ReferenceSpot;
            break;
          }
        }
      }


      /* // Compute the reference profile */
      /* reference_learner_type reference(spec_.const_ref(), data, grid_size_); */


      /* // Get the new columns to set */
      /* af::ref<double> intensity   = data["intensity.prf.value"]; */
      /* af::ref<double> variance    = data["intensity.prf.variance"]; */
      /* af::ref<double> correlation = data["profile_correlation"]; */

      /* // Do the profile fitting for all reflections */
      /* for (std::size_t i = 0; i < data.size(); ++i) { */

      /*   // Initialise some stuff */
      /*   intensity[i] = 0; */
      /*   variance[i] = 0; */
      /*   correlation[i] = 0; */
      /*   flags[i] &= ~IntegratedPrf; */

      /*   // Integrate */
      /*   if (!(flags[i] & DontIntegrate)) { */

      /*     // Get the shoebox */
      /*     const Shoebox<> &sbox = shoebox[i]; */
      /*     DIALS_ASSERT(sbox.is_consistent()); */

      /*     // Get the profile for a given reflection */
      /*     profile_type profile = reference.get( */
      /*         id[i], */
      /*         sbox.panel, */
      /*         s1[i], */
      /*         xyzcal[i][2], */
      /*         sbox.bbox); */

      /*     // Compute the integration mask */
      /*     mask_type mask(sbox.mask.accessor(), false); */
      /*     std::size_t mask_code = Valid | Foreground; */
      /*     std::size_t mask_count = 0; */
      /*     DIALS_ASSERT(profile.accessor().all_eq(sbox.mask.accessor())); */
      /*     for (std::size_t j = 0; j < sbox.mask.size(); ++j) { */
      /*       if ((sbox.mask[j] & mask_code) == mask_code && profile[j] >= 0.0) { */
      /*         mask[j] = true; */
      /*         mask_count++; */
      /*       } */
      /*     } */

      /*     // Perform the profile fit */
      /*     if (mask_count > 0) { */
      /*       fitting_type fit( */
      /*           profile.const_ref(), */
      /*           mask.const_ref(), */
      /*           sbox.data.const_ref(), */
      /*           sbox.background.const_ref()); */

      /*       // Set the data in the reflection */
      /*       intensity[i]   = fit.intensity(); */
      /*       variance[i]    = fit.variance(); */
      /*       correlation[i] = fit.correlation(); */

      /*       // Set the integrated flag */
      /*       flags[i] |= IntegratedPrf; */
      /*     } */
      /*   } */
      /* } */

      // Return the reference learner
      /* return reference; */
    }

  private:

    af::shared<Spec> spec_;
    std::size_t grid_size_;
  };

}} // namespace dials::algorithms


#endif // DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H
