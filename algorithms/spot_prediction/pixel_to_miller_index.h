/*
 * pixel_to_miller_index.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_TO_MILLER_INDEX_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_TO_MILLER_INDEX_H

#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dxtbx::model::BeamBase;
  using dxtbx::model::CrystalBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * A class to transform pixels to miller indices
   */
  class PixelToMillerIndex {
  public:
    /**
     * Initialize with models
     */
    PixelToMillerIndex(const BeamBase &beam,
                       const Detector &detector,
                       const Goniometer &goniometer,
                       const Scan &scan,
                       const CrystalBase &crystal)
        : detector_(detector),
          scan_(scan),
          s0_(beam.get_s0()),
          m2_(goniometer.get_rotation_axis_datum()),
          S_inv_(goniometer.get_setting_rotation().inverse()),
          F_inv_(goniometer.get_fixed_rotation().inverse()),
          A_inv_(crystal.get_A().inverse()) {}

    /**
     * Initialize with models
     */
    PixelToMillerIndex(const BeamBase &beam,
                       const Detector &detector,
                       const CrystalBase &crystal)
        : detector_(detector),
          s0_(beam.get_s0()),
          m2_(0, 0, 0),
          S_inv_(0, 0, 0, 0, 0, 0, 0, 0, 0),
          F_inv_(0, 0, 0, 0, 0, 0, 0, 0, 0),
          A_inv_(crystal.get_A().inverse()) {}

    /**
     * Compute the miller index
     */
    vec3<double> h(std::size_t panel, double x, double y, double z) const {
      DIALS_ASSERT(!(m2_[0] == 0 && m2_[1] == 0 && m2_[2] == 0));

      // Compute the diffracted beam vector
      vec3<double> s1 =
        detector_[panel].get_pixel_lab_coord(vec2<double>(x, y)).normalize()
        * s0_.length();

      // Compute the angle
      double angle = scan_.get_angle_from_array_index(z);

      // Compute the reciprocal lattice vector
      vec3<double> r = s1 - s0_;

      // Create the rotation matrix
      mat3<double> R = scitbx::math::r3_rotation::axis_and_angle_as_matrix(m2_, angle);

      // Compue the miller index
      //  r = S R F A h
      //  where:
      //   S = setting rotation
      //   R = rotation
      //   F = fixed rotation
      //   A = UB
      return A_inv_ * F_inv_ * R.transpose() * S_inv_ * r;
    }

    /**
     * Compute the miller index
     */
    vec3<double> h(std::size_t panel, double x, double y) const {
      // Compute the diffracted beam vector
      vec3<double> s1 =
        detector_[panel].get_pixel_lab_coord(vec2<double>(x, y)).normalize()
        * s0_.length();

      // Compute the reciprocal lattice vector
      vec3<double> r = s1 - s0_;

      // Compue the miller index
      //  r = S R F A h
      //  where:
      //   S = setting rotation
      //   R = rotation
      //   F = fixed rotation
      //   A = UB
      return A_inv_ * r;
    }

    /**
     * Compute the miller index
     */
    vec3<double> q(std::size_t panel, double x, double y) const {
      // Compute the diffracted beam vector
      vec3<double> s1 =
        detector_[panel].get_pixel_lab_coord(vec2<double>(x, y)).normalize()
        * s0_.length();

      // Compute the reciprocal lattice vector
      return s1 - s0_;
    }

  protected:
    Detector detector_;
    Scan scan_;
    vec3<double> s0_;
    vec3<double> m2_;
    mat3<double> S_inv_;
    mat3<double> F_inv_;
    mat3<double> A_inv_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_TO_MILLER_INDEX_H
