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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace shoebox {

  // Use a load of stuff from other namespaces
  using std::floor;
  using std::ceil;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::min;
  using scitbx::af::max;
  using scitbx::af::int6;
  using scitbx::af::double2;
  using scitbx::af::double3;
  using scitbx::af::double4;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  /** Calculate the bounding box for each reflection */
  class BBoxCalculator {

  public:

    BBoxCalculator(const Beam &beam,
                   const Detector &detector,
                   const Goniometer &gonio,
                   const Scan &scan,
                   const af::const_ref<double> &delta_divergence,
                   const af::const_ref<double> &delta_mosaicity)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()),
        detector_(detector),
        scan_(scan),
        delta_divergence_(
          delta_divergence.begin(),
          delta_divergence.end()),
        delta_mosaicity_(
          delta_mosaicity.begin(),
          delta_mosaicity.end()) {
      DIALS_ASSERT(delta_divergence.size() == delta_mosaicity.size());
      DIALS_ASSERT(delta_divergence.size() == scan_.get_num_images());
      DIALS_ASSERT(delta_divergence.all_gt(0.0));
      DIALS_ASSERT(delta_mosaicity.all_gt(0.0));
    }

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
        delta_divergence_(scan.get_num_images(), delta_divergence),
        delta_mosaicity_(scan.get_num_images(), delta_mosaicity) {
      DIALS_ASSERT(delta_divergence > 0.0);
      DIALS_ASSERT(delta_mosaicity > 0.0);
    }

    /**
     * Calculate the bbox on the detector image volume for the reflection.
     *
     * The roi is calculated using the parameters delta_divergence and
     * delta_mosaicity. The reflection mask comprises all pixels where:
     *  |e1| <= delta_d, |e2| <= delta_d, |e3| <= delta_m
     *
     * We transform the coordinates of the box
     *   (-delta_d, -delta_d, 0)
     *   (+delta_d, -delta_d, 0)
     *   (-delta_d, +delta_d, 0)
     *   (+delta_d, +delta_d, 0)
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
      reflection_basis::CoordinateSystem xcs(m2_, s0_, s1, phi);

      // Get the divergence and mosaicity for this point
      int z0 = (int)floor(scan_.get_array_index_from_angle(phi));
      double delta_d = 0, delta_m = 0;
      if (z0 < 0) {
        delta_d = delta_divergence_.front();
        delta_m = delta_mosaicity_.front();
      } else if (z0 >= delta_divergence_.size()) {
        delta_d = delta_divergence_.back();
        delta_m = delta_mosaicity_.back();
      } else {
        delta_d = delta_divergence_[z0];
        delta_m = delta_mosaicity_[z0];
      }

      // Calculate the beam vectors at the following xds coordinates:
      //   (-delta_d, -delta_d, 0)
      //   (+delta_d, -delta_d, 0)
      //   (-delta_d, +delta_d, 0)
      //   (+delta_d, +delta_d, 0)
      double point = delta_d;
      double3 sdash1 = xcs.to_beam_vector(double2(-point, -point));
      double3 sdash2 = xcs.to_beam_vector(double2(+point, -point));
      double3 sdash3 = xcs.to_beam_vector(double2(-point, +point));
      double3 sdash4 = xcs.to_beam_vector(double2(+point, +point));

      // Get the detector coordinates (px) at the ray intersections
      double2 xy1 = detector_[panel].get_ray_intersection_px(sdash1);
      double2 xy2 = detector_[panel].get_ray_intersection_px(sdash2);
      double2 xy3 = detector_[panel].get_ray_intersection_px(sdash3);
      double2 xy4 = detector_[panel].get_ray_intersection_px(sdash4);

      /// Calculate the rotation angles at the following XDS
      // e3 coordinates: -delta_m, +delta_m
      double phi1 = xcs.to_rotation_angle_fast(-delta_m);
      double phi2 = xcs.to_rotation_angle_fast(+delta_m);

      // Get the array indices at the rotation angles
      double z1 = scan_.get_array_index_from_angle(phi1);
      double z2 = scan_.get_array_index_from_angle(phi2);

      // Return the roi in the following form:
      // (minx, maxx, miny, maxy, minz, maxz)
      // Min's are rounded down to the nearest integer, Max's are rounded up
      double4 x(xy1[0], xy2[0], xy3[0], xy4[0]);
      double4 y(xy1[1], xy2[1], xy3[1], xy4[1]);
      double2 z(z1, z2);
      int6 bbox((int)floor(min(x)), (int)ceil(max(x)),
                (int)floor(min(y)), (int)ceil(max(y)),
                (int)floor(min(z)), (int)ceil(max(z)));

      vec2<std::size_t> image_size = detector_[panel].get_image_size();
      std::size_t zrange = delta_divergence_.size();
      bbox[0] = std::max(bbox[0], 0);
      bbox[1] = std::min(bbox[1], (int)image_size[0]);
      bbox[2] = std::max(bbox[2], 0);
      bbox[3] = std::min(bbox[3], (int)image_size[1]);
      bbox[4] = std::max(bbox[4], 0);
      bbox[5] = std::min(bbox[5], (int)zrange);
      return bbox;
    }

    /**
     * Calculate the rois for an array of reflections given by the array of
     * diffracted beam vectors and rotation angles.
     * @param s1 The array of diffracted beam vectors
     * @param phi The array of rotation angles.
     */
    af::shared<int6> operator()(
        const af::const_ref< vec3<double> > &s1,
        const af::const_ref<double> &phi,
        const af::const_ref<std::size_t> &panel) const {
      DIALS_ASSERT(s1.size() == phi.size());
      DIALS_ASSERT(s1.size() == panel.size());
      af::shared<int6> result(s1.size(), af::init_functor_null<int6>());
      for (std::size_t i = 0; i < s1.size(); ++i) {
        result[i] = operator()(s1[i], phi[i], panel[i]);
      }
      return result;
    }

  private:

    vec3<double> s0_;
    vec3<double> m2_;
    Detector detector_;
    Scan scan_;
    af::shared<double> delta_divergence_;
    af::shared<double> delta_mosaicity_;
  };

}}} // namespace dials::algorithms::shoebox

#endif // DIALS_ALGORITHMS_INTEGRATION_BBOX_CALCULATOR_H
