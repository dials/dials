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
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using scitbx::af::mean;
  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * A class to perform profile fitting on reflections.
   */
  class ProfileFittingReciprocalSpace {
  public:

    typedef ReferenceLocator<XdsCircleSampler> locator_type;

    /**
     * Initialise the class. Set the reference profile locator.
     * @param locate The locator function
     */
    ProfileFittingReciprocalSpace(shared_ptr<locator_type> locate)
     : locate_(locate) {}

    /**
     * Perform the profile fitting on all the reflections
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
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
      flex_double c = reflection.get_transformed_shoebox();

      // Get the transformed background
      // HACK ALERT! Should fix: setting to mean of shoebox background
      flex_double b = flex_double(c.accessor(),
        mean(reflection.get_shoebox_background().const_ref()));

      // Get the reference profile at the reflection coordinate
      double3 coord(reflection.get_image_coord_px()[0],
                    reflection.get_image_coord_px()[1],
                    reflection.get_frame_number());
      flex_double p = locate_->profile(coord);

      // Do the profile fitting and set the intensity and variance
      ProfileFitting2 fit(p, c, b);
      reflection.set_intensity(fit.intensity());
      reflection.set_intensity_variance(fit.variance());
    }

  private:
    shared_ptr<locator_type> locate_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_RECIPROCAL_SPACE_H */
