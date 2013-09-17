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
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using scitbx::af::mean;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dials::model::Reflection;

  /**
   * A class to do 3D summation integration in reciprocal space
   */
  class SummationReciprocalSpace {
  public:

    typedef Summation integrator;

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
      af::const_ref< double, af::c_grid<3> > c =
        r.get_transformed_shoebox().const_ref();

      // Get the transformed background
      // HACK ALERT! Should fix: setting to mean of shoebox background
      af::versa< double, af::c_grid<3> > b(c.accessor(), mean(
        r.get_shoebox_background().const_ref()));

      // Integrate the reflection
      integrator result = integrator(c, b.const_ref());

      // Set intensity data in reflection container
      r.set_intensity(result.intensity());
      r.set_intensity_variance(result.variance());
    }

    /**
     * Integrate a list of reflections
     * @param reflections The reflection list
     */
    void operator()(af::ref<Reflection> reflections) const {
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
