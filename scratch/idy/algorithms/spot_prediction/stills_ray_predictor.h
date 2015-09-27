/*
 * stills_ray_predictor.h
 *
 *  Copyright (C) 2015 Diamond Light Source, CCP4
 *
 *  Author: Nicholas Sauter
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef IDY_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
#define IDY_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H

#include <cmath>
#include <cctbx/miller.h>
#include <scitbx/array_family/small.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/ray.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using std::sqrt;
  using scitbx::vec3;
  using scitbx::mat3;
  using dials::model::Ray;

  /**
   *  Parent class StillsRayPredictor always assumes that the center of impact is at
   *  radius r from the origin of reciprocal space
   *
   *  In contrast StillsCentralImpactRayPredictor takes the picture from
   *  Juers Fig. 4 (2007) Acta D 63, 1139-1153.  If the angular separation
   *  between Ewald sphere surface and relp is greater than omega/2 (as
   *  defined in that Figure) then take the ray to [d/s + da/a] sphere
   *  drawn at omega/2.
   */
  class StillsCentralImpactRayPredictor {
  public:

    StillsCentralImpactRayPredictor(vec3<double> s0, double const& mos)
      : s0_(s0), half_mosaicity_rad_(mos) {
      SCITBX_ASSERT (half_mosaicity_rad_>=0.);
      DIALS_ASSERT(s0_.length() > 0.0);
      unit_s0_ = s0_.normalize();
    }

    Ray operator()(miller_index h, mat3<double> ub) {

      // Calculate the reciprocal space vector and required unit vectors
      vec3<double> q = ub * h;
      vec3<double> e1 = q.cross(unit_s0_).normalize();
      vec3<double> c0 = unit_s0_.cross(e1).normalize();

      // Calculate the vector rotated to the Ewald sphere
      double qq = q.length_sq();
      double lambda = 1. / s0_.length();
      double a = 0.5 * qq * lambda;
      double tmp = qq - a*a;
      DIALS_ASSERT(tmp > 0.0);
      double b = std::sqrt(tmp);
      vec3<double> r = -1.0 * a * unit_s0_ + b * c0;

      // Calculate delpsi value
      vec3<double> q0 = q.normalize();
      vec3<double> q1 = q0.cross(e1).normalize();
      delpsi_ = -1.0 * atan2(r*q1, r*q0);

      // Calculate the Ray (default zero angle and 'entering' as false)
      vec3<double> s1 = (s0_ + r).normalize() * s0_.length();
      if (std::abs(delpsi_) > half_mosaicity_rad_) {
        // Use rodrigues rotation formula to rotate Q about e1.
        vec3<double> r_new = q.unit_rotate_around_origin(
        //xxxxxxx does not seem to work or have the desired effect
                             e1, (delpsi_>0.)?10*half_mosaicity_rad_:-10*half_mosaicity_rad_);
        s1 = (s0_ + r_new).normalize() * s0_.length();
      }
      return Ray(s1, 0.0, false);
    }
    double get_delpsi() const {
      return delpsi_;
    }

  protected:
    vec3<double> s0_;
    vec3<double> unit_s0_;
    double delpsi_;
    double half_mosaicity_rad_;

  };

}} // namespace dials::algorithms

#endif // IDY_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
