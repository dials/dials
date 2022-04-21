/*
 * stills_ray_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source, CCP4
 *
 *  Author: David Waterman
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H

#include <cmath>
#include <cctbx/miller.h>
#include <scitbx/array_family/small.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/ray.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::model::Ray;
  using scitbx::mat3;
  using scitbx::vec3;
  using std::sqrt;

  /**
   *  Predict for a relp based on the current states of models of the
   *  experimental geometry. Here we assume the crystal UB already puts hkl in
   *  reflecting position, so no rotation is required.
   *
   *  Generally, this assumption is not true: most relps are not touching the
   *  Ewald sphere on a still image, but are observed due to their finite
   *  mosaicity and the finite width of the Ewald sphere. Rather than employing
   *  a complicated polychromatic and mosaic model, here we take a naive
   *  approach, which is likely to introduce (small?) errors in the direction of
   *  predicted rays.
   */
  class StillsRayPredictor {
  public:
    typedef cctbx::miller::index<> miller_index;

    StillsRayPredictor(vec3<double> s0) : s0_(s0) {
      DIALS_ASSERT(s0_.length() > 0.0);
      unit_s0_ = s0_.normalize();
    }

    Ray operator()(miller_index h, mat3<double> ub) {
      // Calculate the reciprocal space vector and required unit vectors
      vec3<double> q = ub * h;
      DIALS_ASSERT(q.length() > 0);
      vec3<double> e1 = q.cross(unit_s0_).normalize();
      vec3<double> c0 = unit_s0_.cross(e1).normalize();

      // Calculate the vector rotated to the Ewald sphere
      double qq = q.length_sq();
      double lambda = 1. / s0_.length();
      double a = 0.5 * qq * lambda;
      double tmp = qq - a * a;
      DIALS_ASSERT(tmp > 0.0);
      double b = std::sqrt(tmp);
      vec3<double> r = -1.0 * a * unit_s0_ + b * c0;

      // Calculate delpsi value
      vec3<double> q0 = q.normalize();
      vec3<double> q1 = q0.cross(e1).normalize();
      delpsi_ = -1.0 * atan2(r * q1, r * q0);

      // Calculate the Ray (default zero angle and 'entering' as false)
      vec3<double> s1 = (s0_ + r).normalize() * s0_.length();
      return Ray(s1, 0.0, false);
    }

    double get_delpsi() const {
      return delpsi_;
    }

  private:
    vec3<double> s0_;
    vec3<double> unit_s0_;
    double delpsi_;
  };

  /**
   * Alternate version of StillsRayPredictor in which no rotation to the
   * surface of the Ewald sphere is assumed to take place. Rather, the relp
   * is assumed to have spherical extent and the ray is directed towards the
   * centre of the region of intersection between the relp sphere and the Ewald
   * sphere.
   *
   * Currently, DeltaPsi is still calculated and will be used as a restraint
   * in refinement in the same way as for StillsRayPredictor. This permits
   * easy comparison between the two versions. However, it may be more
   * appropriate to use the distance between the centre of the relp and the
   * Ewald sphere, epsilon, as a restraint in future.
   */
  class SphericalRelpStillsRayPredictor {
  public:
    typedef cctbx::miller::index<> miller_index;

    SphericalRelpStillsRayPredictor(vec3<double> s0) : s0_(s0) {
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
      double s0len = s0_.length();
      double a = 0.5 * qq / s0len;
      double tmp = qq - a * a;
      DIALS_ASSERT(tmp > 0.0);
      double b = std::sqrt(tmp);
      vec3<double> r = -1.0 * a * unit_s0_ + b * c0;

      // Calculate delpsi value
      vec3<double> q0 = q.normalize();
      vec3<double> q1 = q0.cross(e1).normalize();
      delpsi_ = -1.0 * atan2(r * q1, r * q0);

      // Calculate epsilon value
      double qs0 = q * s0_;
      tmp = qq + 2.0 * qs0 + s0len * s0len;
      DIALS_ASSERT(tmp > 0.0);
      epsilon_ = std::sqrt(tmp) - s0len;

      // Calculate the Ray (default zero angle and 'entering' as false)
      vec3<double> s1 = (s0_ + q).normalize() * s0_.length();
      return Ray(s1, 0.0, false);
    }

    double get_delpsi() const {
      return delpsi_;
    }

    double get_epsilon() const {
      return epsilon_;
    }

  private:
    vec3<double> s0_;
    vec3<double> unit_s0_;
    double delpsi_;
    double epsilon_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
