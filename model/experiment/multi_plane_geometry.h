/*
 * multi_plane_geometry.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_MULTI_PLANE_GEOMETRY_H
#define DIALS_MODEL_EXPERIMENT_MULTI_PLANE_GEOMETRY_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double4;

  // Create flex array typedefs
  typedef scitbx::af::flex<mat3<double> >::type flex_mat3_double;
  typedef scitbx::af::flex<double4>::type flex_double4;

  /**
   * A functor to calculate the intersection of a diffracted beam vector coming
   * from the origin of the laboratory coordinate system with multiple planes. 
   * The planes are given by the D matrix which is the inverse of the matrix
   * 
   * d = [[dx1, dy1, d01], [dx2, dy2, d02], [dx3, dy3, d03]]
   *
   * where the dx, dy, d0 are the basis vectors of the plane.
   *
   * The coordinate returned is in the units given by the basis vectors.
   *
   * @todo This functor currently uses a brute force approach to finding the
   * detector panel which contains the reflection. All panels are checked
   * and where more than one panel could record the reflection, the closest
   * is chosen. This is ok for a small number of panels but doesn't scale
   * well for a large number.
   */
  class BeamMultiPlaneIntersection {
  public:
    typedef std::pair<int, vec2<double> > coord_type;
  
    /**
     * Initialise the transform from the D matrices and plane extents
     * @param D The D matrix array
     * @param extents The extents of the planes
     */
    BeamMultiPlaneIntersection(const flex_mat3_double &D, 
        const flex_double4 &extents)
      : D_(D),
        extents_(extents) {}

    /**
     * Calculate the intersection point on the planes.
     * @param s1 The diffracted beam vector
     * @returns The (panel, (x, y)) detector coordinate
     */
    coord_type operator()(vec3 <double> s1) const {
      coord_type pxy(-1, vec2<double>(0, 0));
      double w_max = 0;

      // Loop through all detectors. If the w component of the (u, v, w)
      // vector points in the correct direction, then calculate the coordinate.
      // If the coordinate is valid and the w component is greater than that of
      // the current closest valid coordinate, then set this coordinate as the
      // current best bet.
      for (std::size_t i = 0; i < D_.size(); ++i) {
        vec3 <double> v = D_[i] * s1;
        if (v[2] > 0) {
          coord_type pxy_temp(i, vec2<double>(v[0] / v[2], v[1] / v[2]));
          if (is_coord_valid(pxy_temp) && v[2] > w_max) {
            pxy = pxy_temp;
            w_max = v[2];
          }
        }
      }

      // If no coordinate was found then raise an exception
      // otherwise return the coordinate.
      DIALS_ASSERT(w_max > 0);
      return pxy;
    }

    /**
     * Calculate the intersection point on a given detector plane.
     * @param s1 The diffracted beam vector
     * @param plane The plane number
     * @returns The (x, y) detector coordinate
     */
    vec2 <double> operator()(vec3 <double> s1, int plane) const {
      vec3 <double> v = D_[plane] * s1;
      DIALS_ASSERT(v[2] > 0);
      return vec2 <double>(v[0] / v[2], v[1] / v[2]);
    }

  private:
    /** Check if a coordinate is valid */
    bool is_coord_valid(int plane, vec2<double> xy) const {
      return extents_[plane][0] <= xy[0] && xy[0] < extents_[plane][2]
          && extents_[plane][1] <= xy[1] && xy[1] < extents_[plane][3];
    }

    /** Check if a coordinate is valid */
    bool is_coord_valid(coord_type pxy) const {
      return is_coord_valid(pxy.first, pxy.second);
    }

    flex_mat3_double D_;
    flex_double4 extents_;
  };
  
  /**
   * A functor to calculate the laboratory coordinate of a point on the plane.
   * The plane is specified by the d matrix
   * 
   * d = [[dx1, dy1, d01], [dx2, dy2, d02], [dx3, dy3, d03]]
   *
   * where the dx, dy, d0 are the basis vectors of the plane.
   */
  class MultiPlaneToLabTransform {
  public:
    typedef std::pair<int, vec2<double> > coord_type;
    
    /**
     * Initialise the transform with the d matrix array
     * @param d The d matrix
     */
    MultiPlaneToLabTransform(const flex_mat3_double &d)
      : d_(d) {}
      
    /**
     * Calculate the lab coordinate of the point on the plane
     * @param plane The plane number
     * @param xy The plane coordinate
     * @returns The (x, y, z) laboratory coordinate
     */
    vec3 <double> operator()(int plane, vec2<double> xy) const {
      vec3<double> v(xy[0], xy[1], 1.0);
      return d_[plane] * v;
    }

    /**
     * Calculate the lab coordinate of the point on the plane
     * @param pxy The plane number and plane coordinate
     * @returns The (x, y, z) laboratory coordinate
     */    
    vec3<double> operator()(coord_type pxy) const {
      return operator()(pxy.first, pxy.second);
    }
    
  private:
    flex_mat3_double d_;  
  };  
}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_MULTI_PLANE_GEOMETRY_H
