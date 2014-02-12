/*
 * stills_ray_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
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

namespace dials { namespace algorithms {

  using scitbx::vec3;
  using scitbx::mat3;

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
      : s0_(s0) {}

    af::small<Ray,2> operator()(miller_index h, mat3<double> ub) const {

      // Calculate the reciprocal space vector
      vec3<double> r = ub * h;
      vec3<double> s1 = (s0_ + r).normalize() * s0_.length();

      // Add both enter and exit reflections
      af::small<Ray,2> result;
      result.push_back(Ray(s1, 0.0, true));
      result.push_back(Ray(s1, 0.0, false));
      return result;
    }

  private:
    vec3<double> s0_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_STILLS_RAY_PREDICTOR_H
