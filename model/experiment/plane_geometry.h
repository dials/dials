/*
 * plane_geometry.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_PLANE_GEOMETRY_H
#define DIALS_MODEL_EXPERIMENT_PLANE_GEOMETRY_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;

  /**
   * A functor to calculate the intersection of a diffracted beam vector coming
   * from the origin of the laboratory coordinate system with a plane. The
   * plane is given by the D matrix which is the inverse of the matrix
   * 
   * d = [[dx1, dy1, d01], [dx2, dy2, d02], [dx3, dy3, d03]]
   *
   * where the dx, dy, d0 are the basis vectors of the plane.
   *
   * The coordinate returned is in the units given by the basis vectors.
   */
  class BeamPlaneIntersection {
  public:
    /**
     * Initialise the transform with the D matrix
     * @param D The D matrix
     */
    BeamPlaneIntersection(mat3 <double> D)
      : D_(D) {}

    /**
     * Calculate the intersection point on the plane.
     * @param s1 The diffracted beam vector
     * @returns The (x, y) detector coordinate
     */
    vec2 <double> operator()(vec3 <double> s1) const {
      vec3 <double> v = D_ * s1;
      DIALS_ASSERT(v[2] > 0);
      return vec2<double>(v[0] / v[2], v[1] / v[2]);
    }
        
  private:
    mat3 <double> D_;
  };

  /**
   * A functor to calculate the laboratory coordinate of a point on the plane.
   * The plane is specified by the d matrix
   * 
   * d = [[dx1, dy1, d01], [dx2, dy2, d02], [dx3, dy3, d03]]
   *
   * where the dx, dy, d0 are the basis vectors of the plane.
   */
  class PlaneToLabTransform {
  public:
    /**
     * Initialise the transform with the d matrix
     * @param d The d matrix
     */
    PlaneToLabTransform(mat3<double> d)
      : d_(d) {}
      
    /**
     * Calculate the lab coordinate of the point on the plane
     * @param xy The plane coordinate
     * @returns The (x, y, z) laboratory coordinate
     */
    vec3 <double> operator()(vec2<double> xy) const {
      vec3<double> v(xy[0], xy[1], 1.0);
      return d_ * v;
    }
    
  private:
    mat3 <double> d_;  
  };

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_PLANE_GEOMETRY_H
