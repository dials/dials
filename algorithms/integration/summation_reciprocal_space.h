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

      // Integrate the reflection
      integrator result = integrator(r.get_transformed_shoebox());

      // Set intensity data in reflection container
      r.set_intensity(result.intensity());
      r.set_intensity_variance(result.variance());

      // Get the centroid information from the result
      vec3<double> centroid = result.centroid();
      vec3<double> variance = result.centroid_standard_error_sq();
      vec3<double> sq_width = result.centroid_variance();

      // Get stuff from reflection struct
      int panel = r.get_panel_number();
      vec3<double> s1 = r.get_beam_vector();
      double phi = r.get_rotation_angle();

      // Initialise the transform from reciprocal space
      XdsCoordinateSystem xcs(s0_, s1, m2_, phi);
      FromXdsToBeamVector transform_s1(xcs, s1);
      FromXdsE3ToPhi transform_phi(xcs.get_zeta(), phi);

      // Get the mm centroid from the reciprocal space centroid
      vec3<double> s1_centroid = transform_s1(centroid);
      vec2<double> mm_centroid = (*detector_)[panel].get_ray_intersection(
        s1_centroid);

      // Get the phi centroid from the reciprocal space centroid
      double phi_centroid = transform_phi(centroid[2]);

      // Set the centroid data in the reflection container
      r.set_centroid_position(vec3<double>(
        mm_centroid[0], mm_centroid[1], phi_centroid));
    }

    /**
     * Integrate a list of reflections
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          this->operator()(reflections[i]);
        } catch (dials::error) {
          continue;
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
