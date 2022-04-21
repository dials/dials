/*
 * ray.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_RAY_H
#define DIALS_MODEL_DATA_RAY_H

#include <scitbx/vec3.h>

namespace dials { namespace model {

  using scitbx::vec3;

  /**
   * A class to represent a ray.
   */
  struct Ray {
    vec3<double> s1;
    double angle;
    bool entering;

    Ray() {}

    Ray(vec3<double> s1_, double angle_, bool entering_)
        : s1(s1_), angle(angle_), entering(entering_) {}
  };

}}  // namespace dials::model

#endif  // DIALS_MODEL_DATA_RAY_H
