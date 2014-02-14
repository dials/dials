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
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/reflection_basis/transform.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using dials::model::Shoebox;
  using dials::algorithms::reflection_basis::transform::Forward;

  /**
   * A class to perform profile fitting on reflections.
   */
  class ProfileFittingReciprocalSpace {
  public:

    typedef Shoebox<>::float_type FloatType;
    typedef ReferenceLocator<FloatType, XdsCircleSampler> locator_type;

    /**
     * Initialise the class. Set the reference profile locator.
     * @param locate The locator function
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
    af::shared< vec2<double> > operator()(
        const af::const_ref< Forward<> > &profiles,
        const af::const_ref< vec3<double> > &coords) const {
      DIALS_ASSERT(profiles.size() == coords.size());
      af::shared< vec2<double> > result(profiles.size());
      for (std::size_t i = 0; i < profiles.size(); ++i) {
        try {
          result[i] = this->operator()(profiles[i], coords[i]);
        } catch (dials::error) {
          result[i] = vec2<double>(0.0, -1.0);
        }
      }
      return result;
    }

    /**
     * Perform the profile fitting on a reflection
     * @param reflection The reflection to process
     */
    vec2<double> operator()(const Forward<> &profile, vec3<double> coord) const {

      typedef af::versa < FloatType, af::c_grid<3> > profile_type;
      typedef af::const_ref< FloatType, af::c_grid<3> > profile_ref_type;

      // Get the transformed shoebox
      profile_ref_type c = profile.profile().const_ref();
      profile_ref_type b = profile.background().const_ref();
      profile_type p = locate_->profile(coord);

      // Do the profile fitting and set the intensity and variance
      ProfileFitting<FloatType> fit(p.const_ref(), c, b, max_iter_);
      DIALS_ASSERT(fit.niter() < max_iter_);
      return vec2<double>(fit.intensity(), fit.variance());
    }

  private:
    shared_ptr<locator_type> locate_;
    std::size_t max_iter_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H
