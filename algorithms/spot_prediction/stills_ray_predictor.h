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

#include <cctbx/miller.h>
#include <scitbx/array_family/small.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/ray.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;
  using scitbx::mat3;
  using dials::model::Ray;

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

    StillsRayPredictor(vec3<double> s0)
      : s0_(s0) {
      DIALS_ASSERT(s0_.length() > 0.0);
    }

    Ray operator()(miller_index h, mat3<double> ub) {

      // Calculate the reciprocal space vector
      vec3<double> r = ub * h;
      vec3<double> s1 = (s0_ + r).normalize() * s0_.length();

      // Calculate delpsi value
      delpsi_ = 0.0; // set dummy value for now

      // Calculate the Ray (default zero angle and 'entering' as true)
      return Ray(s1, 0.0, true);
    }

    double get_delpsi() const {
      return delpsi_;
    }

  private:
    vec3<double> s0_;
    double delpsi_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
