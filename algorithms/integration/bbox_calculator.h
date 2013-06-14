/*
 * bbox_calculator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_BBOX_CALCULATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_BBOX_CALCULATOR_H

#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/model/data/reflection.h>
#include "xds_coordinate_system.h"
#include "from_xds_to_beam_vector.h"
#include "from_xds_e3_to_phi.h"

namespace dials { namespace algorithms {

  // Use a load of stuff from other namespaces
  using std::floor;
  using std::ceil;
  using scitbx::vec3;
  using scitbx::af::min;
  using scitbx::af::max;
  using scitbx::af::int6;
  using scitbx::af::double2;
  using scitbx::af::double3;
  using scitbx::af::double4;
  using scitbx::af::flex_double;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;
  typedef scitbx::af::flex<int6>::type flex_int6;

  /** Calculate the bounding box for each reflection */
  class BBoxCalculator {

  public:

    /**
     * Initialise the bounding box calculation.
     * @param beam The beam parameters
     * @param detector The detector parameters
     * @param goniometer The goniometer parameters
     * @param delta_divergence The xds delta_divergence parameter
     * @param delta_mosaicity The xds delta_mosaicity parameter
     */
    BBoxCalculator(const Beam &beam,
                   const Detector &detector,
                   const Goniometer &gonio,
                   const Scan &scan,
                   double delta_divergence,
                   double delta_mosaicity)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()),
        detector_(detector),
        scan_(scan),
        delta_divergence_(delta_divergence),
        delta_mosaicity_(delta_mosaicity) {
      DIALS_ASSERT(delta_divergence > 0.0);
      DIALS_ASSERT(delta_mosaicity > 0.0);
    }

    /**
     * Calculate the bbox on the detector image volume for the reflection.
     *
     * The roi is calculated using the parameters delta_divergence and
     * delta_mosaicity. The reflection mask comprises all pixels where:
     *  |e1| <= delta_d / 2, |e2| <= delta_d < 2, |e3| <= delta_m / 2
     *
     * We transform the coordinates of the box
     *   (-delta_d/2, -delta_d/2, 0)
     *   (+delta_d/2, -delta_d/2, 0)
     *   (-delta_d/2, +delta_d/2, 0)
     *   (+delta_d/2, +delta_d/2, 0)
     *
     * to the detector image volume and return the minimum and maximum values
     * for the x, y, z image volume coordinates.
     *
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     * @returns A 6 element array: (minx, maxx, miny, maxy, minz, maxz)
     */
    int6 operator()(vec3 <double> s1, double phi, std::size_t panel) const {

      // Ensure our values are ok
      DIALS_ASSERT(s1.length_sq() > 0);

      // Create the coordinate system for the reflection
      XdsCoordinateSystem xcs(s0_, s1, m2_, phi);

      // Create the transformer from the xds coordinate system to the detector
      FromXdsToBeamVector calculate_beam_vector(xcs, s1);

      // Create the transformer from Xds E3 to rotation angle
      FromXdsE3ToPhi calculate_rotation_angle(xcs.get_zeta(), phi);

      // Calculate the beam vectors at the following xds coordinates:
      //   (-delta_d/2, -delta_d/2, 0)
      //   (+delta_d/2, -delta_d/2, 0)
      //   (-delta_d/2, +delta_d/2, 0)
      //   (+delta_d/2, +delta_d/2, 0)
      double point = delta_divergence_ / 2.0;
      double3 sdash1 = calculate_beam_vector(double3(-point, -point, 0.0));
      double3 sdash2 = calculate_beam_vector(double3(+point, -point, 0.0));
      double3 sdash3 = calculate_beam_vector(double3(-point, +point, 0.0));
      double3 sdash4 = calculate_beam_vector(double3(+point, +point, 0.0));

      // Get the detector coordinates (px) at the ray intersections
      double2 xy1 = detector_[panel].get_ray_intersection_px(sdash1);
      double2 xy2 = detector_[panel].get_ray_intersection_px(sdash2);
      double2 xy3 = detector_[panel].get_ray_intersection_px(sdash3);
      double2 xy4 = detector_[panel].get_ray_intersection_px(sdash4);

      /// Calculate the rotation angles at the following XDS
      // e3 coordinates: -delta_m/2, +delta_m/2
      double phi1 = calculate_rotation_angle(-delta_mosaicity_ / 2.0);
      double phi2 = calculate_rotation_angle(+delta_mosaicity_ / 2.0);

      // Get the array indices at the rotation angles
      double z1 = scan_.get_array_index_from_angle(phi1);
      double z2 = scan_.get_array_index_from_angle(phi2);

      // Return the roi in the following form:
      // (minx, maxx, miny, maxy, minz, maxz)
      // Min's are rounded down to the nearest integer, Max's are rounded up
      double4 x(xy1[0], xy2[0], xy3[0], xy4[0]);
      double4 y(xy1[1], xy2[1], xy3[1], xy4[1]);
      double2 z(z1, z2);
      return int6 (
        (int)floor(min(x)), (int)ceil(max(x)),
        (int)floor(min(y)), (int)ceil(max(y)),
        (int)floor(min(z)), (int)ceil(max(z)));
    }

    /**
     * Calculate the rois for an array of reflections given by the array of
     * diffracted beam vectors and rotation angles.
     * @param s1 The array of diffracted beam vectors
     * @param phi The array of rotation angles.
     */
    flex_int6 operator()(const flex_vec3_double &s1,
        const flex_double &phi, std::size_t panel) const {
      DIALS_ASSERT(s1.size() == phi.size());
      flex_int6 result(s1.size());
      for (std::size_t i = 0; i < s1.size(); ++i) {
        result[i] = operator()(s1[i], phi[i], panel);
      }
      return result;
    }

    /**
     * Calculate the rois for a reflection given by reflection struct
     * @param reflection The reflection data
     */
    void operator()(Reflection &reflection) const {
      reflection.set_bounding_box(
        operator()(reflection.get_beam_vector(),
                   reflection.get_rotation_angle(),
                   reflection.get_panel_number()));
    }

    /**
     * Calculate the rois for an array of reflections given by the array of
     * reflections
     * @param reflections The list of reflections
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        operator()(reflections[i]);
      }
    }

  private:

    vec3<double> s0_;
    vec3<double> m2_;
    Detector detector_;
    Scan scan_;
    double delta_divergence_;
    double delta_mosaicity_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_BBOX_CALCULATOR_H
