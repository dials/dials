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

#ifndef DIALS_ALGORITHMS_INTEGRATION_FIT_FIT_H
#define DIALS_ALGORITHMS_INTEGRATION_FIT_FIT_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/integration/profile/fitting.h>


namespace dials { namespace algorithms {

  using af::DontIntegrate;
  using af::IntegratedPrf;
  using model::Shoebox;
  using model::Valid;
  using model::Foreground;


  /**
   * A class to compute a single reference profile.
   */
  template <typename FloatType>
  class SingleReferenceProfile {
  public:

    typedef af::versa< FloatType, af::c_grid<3> > profile_type;

    SingleReferenceProfile(af::reflection_table data) {

    }

    profile_type get(
        std::size_t id,
        std::size_t panel,
        vec3<double> xyz,
        int6 bbox) {
      return profile_type();
    }

  };


  /**
   * A class to perform image space profile fitting.
   */
  class ImageSpaceProfileFitting {
  public:

    typedef Shoebox<>::float_type float_type;
    typedef SingleReferenceProfile<float_type> reference_type;
    typedef reference_type::profile_type profile_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;
    typedef ProfileFitting<float_type> fitting_type;

    ImageSpaceProfileFitting() {}

    void exec(af::reflection_table data) const {

      // Check the input contains expected fields
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("xyzcal.px"));
      DIALS_ASSERT(data.contains("flags"));

      // Compute the reference profile
      reference_type reference(data);

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
        if (!(flags[i] & DontIntegrate)) {

          // Get the shoebox
          const Shoebox<> &sbox = shoebox[i];

          // Get the profile for a given reflection
          profile_type profile = reference.get(
              id[i],
              sbox.panel,
              xyzcal[i],
              sbox.bbox);

          // Compute the integration mask
          mask_type mask(sbox.mask.accessor(), false);
          std::size_t mask_code = Valid | Foreground;
          for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
            mask[j] = (sbox.mask[j] & mask_code) == mask_code;
          }

          // Perform the profile fit
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

        } else {
          intensity[i] = 0;
          variance[i] = 0;
          correlation[i] = 0;
          flags[i] &= ~IntegratedPrf;
        }
      }

    }
  };

}} // namespace dials::algorithms


#endif // DIALS_ALGORITHMS_INTEGRATION_FIT_FIT_H
