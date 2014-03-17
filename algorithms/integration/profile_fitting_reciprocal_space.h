/*
 * profile_fitting_reciprocal_space.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H

#include <omptbx/omp_or_stubs.h>
#include <boost/shared_ptr.hpp>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/algorithms/integration/profile/fitting.h>
#include <dials/algorithms/integration/profile/grid_sampler.h>
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/transformed_shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using dials::model::Shoebox;
  using dials::model::TransformedShoebox;

  /**
   * A class to perform profile fitting on reflections.
   */
  class ProfileFittingReciprocalSpace {
  public:

    typedef Shoebox<>::float_type FloatType;
    typedef ReferenceLocator<FloatType, GridSampler> locator_type;

    /**
     * Initialise the class. Set the reference profile locator.
     * @param locate The locator function
     * @param mask_range n_sigma_mask / n_sigma_grid
     */
    ProfileFittingReciprocalSpace(
        shared_ptr<locator_type> locate,
        std::size_t max_iter)
     : locate_(locate),
       max_iter_(max_iter) {
      DIALS_ASSERT(locate != NULL);
      DIALS_ASSERT(max_iter_ > 0);
    }

    /**
     * Perform the profile fitting on all the reflections
     * @param reflections The reflection list
     */
    af::shared< vec3<double> > operator()(
        const af::const_ref<TransformedShoebox> &profiles,
        const af::const_ref< vec3<double> > &coords) const {
      DIALS_ASSERT(profiles.size() == coords.size());
      af::shared< vec3<double> > result(profiles.size());
      for (std::size_t i = 0; i < profiles.size(); ++i) {
        try {
          result[i] = this->operator()(profiles[i], coords[i]);
        } catch (dials::error) {
          result[i] = vec3<double>(0.0, -1.0, 0.0);
        }
      }
      return result;
    }

    /**
     * Perform the profile fitting on a reflection
     * @param reflection The reflection to process
     */
    vec3<double> operator()(const TransformedShoebox &profile, vec3<double> coord) const {

      typedef af::versa < FloatType, af::c_grid<3> > profile_type;
      typedef af::versa < bool, af::c_grid<3> > profile_mask_type;
      typedef af::const_ref< FloatType, af::c_grid<3> > profile_ref_type;

      // Get the transformed shoebox
      profile_ref_type c = profile.data.const_ref();
      profile_ref_type b = profile.background.const_ref();
      profile_type p = locate_->profile(coord);
      profile_mask_type m = locate_->mask(coord);

      // Do the profile fitting and set the intensity and variance
      ProfileFitting<FloatType> fit(p.const_ref(), m.const_ref(), c, b, 1e-3, max_iter_);
      DIALS_ASSERT(fit.niter() < max_iter_);
      return vec3<double>(fit.intensity(), fit.variance(), fit.correlation());
    }

  private:
    shared_ptr<locator_type> locate_;
    std::size_t max_iter_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H
