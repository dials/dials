/*
 * summation_reciprocal_space.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_SUMMATION_RECIPROCAL_SPACE_H
#define DIALS_ALGORITHMS_INTEGRATION_SUMMATION_RECIPROCAL_SPACE_H

#include <boost/shared_ptr.hpp>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dials/algorithms/integration/summation.h>
#include <dials/algorithms/integration/xds_coordinate_system.h>
#include <dials/algorithms/integration/from_xds_to_beam_vector.h>
#include <dials/algorithms/integration/from_xds_e3_to_phi.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using scitbx::af::mean;
  using scitbx::af::flex_double;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * A class to do 3D summation integration in reciprocal space
   */
  class SummationReciprocalSpace {
  public:

    typedef IntegrateBySummation integrator;

    /** Init the algorithm. */
    SummationReciprocalSpace(const shared_ptr<Beam> &beam,
                             const shared_ptr<Detector> &detector,
                             const shared_ptr<Goniometer> &gonio)
      : s0_(beam->get_s0()),
        m2_(gonio->get_rotation_axis()),
        detector_(detector) {}

    /**
     * Integrate a reflection
     * @param r The reflection container
     */
    void operator()(Reflection &r) const {

      // Get the transformed shoebox
      flex_double c = r.get_transformed_shoebox();

      // Get the transformed background
      // HACK ALERT! Should fix: setting to mean of shoebox background
      flex_double b = flex_double(c.accessor(),
        mean(r.get_shoebox_background().const_ref()));

      // Integrate the reflection
      integrator result = integrator(c, b);

      // Set intensity data in reflection container
      r.set_intensity(result.intensity());
      r.set_intensity_variance(result.variance());
    }

    /**
     * Integrate a list of reflections
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          if (reflections[i].is_valid()) {
            this->operator()(reflections[i]);
          }
        } catch (dials::error) {
          reflections[i].set_valid(false);
        }
      }
    }

  private:
    vec3<double> s0_;
    vec3<double> m2_;
    shared_ptr<Detector> detector_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_SUMMATION_RECIPROCAL_SPACE_H */
