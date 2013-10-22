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
#include <dials/model/data/reflection.h>
#include <dials/algorithms/integration/profile/fitting.h>
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using dials::model::Reflection;

  /**
   * A class to perform profile fitting on reflections.
   */
  class ProfileFittingReciprocalSpace {
  public:

    typedef Reflection::float_type FloatType;
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
      DIALS_ASSERT(max_iter_ > 0);
    }

    /**
     * Perform the profile fitting on all the reflections
     * @param reflections The reflection list
     */
    void operator()(af::ref<Reflection> reflections) const {
      #pragma omp parallel for
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].is_valid()) {
          try {
            this->operator()(reflections[i]);
          } catch (dials::error) {
            reflections[i].set_valid(false);
          }
        }
      }
    }

    /**
     * Perform the profile fitting on a reflection
     * @param reflection The reflection to process
     */
    void operator()(Reflection &reflection) const {

      // Get the transformed shoebox
      af::const_ref<FloatType, af::c_grid<3> > c =
        reflection.get_transformed_shoebox().const_ref();
      af::const_ref<FloatType, af::c_grid<3> > b =
        reflection.get_transformed_shoebox_background().const_ref();

      // Get the reference profile at the reflection coordinate
      double3 coord(reflection.get_image_coord_px()[0],
                    reflection.get_image_coord_px()[1],
                    reflection.get_frame_number());
      af::versa<FloatType, af::c_grid<3> > p = locate_->profile(coord);

      // Do the profile fitting and set the intensity and variance
      ProfileFitting<FloatType> fit(p.const_ref(), c, b, max_iter_);
      DIALS_ASSERT(fit.niter() < max_iter_);
      reflection.set_intensity(fit.intensity());
      reflection.set_intensity_variance(fit.variance());
    }

  private:
    shared_ptr<locator_type> locate_;
    std::size_t max_iter_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H
