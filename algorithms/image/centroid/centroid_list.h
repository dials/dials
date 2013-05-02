/*
 * centroid_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H
#define DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::flex_int;

  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;
  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;

  template <typename T>
  T sqr(const T &a) {
    return a * a;
  }

  /**
   * Class to calculate the 2D centroid of a list of coordinates
   */
  class CentroidList2d {
  public:

    /**
     * Calculate the centroid.
     * @param coord The list of coordinates.
     * @param value The list of values
     */
    CentroidList2d(const flex_int &value, const flex_vec2_double &coord)
      : counts_(0),
        position_(0, 0),
        variance_(0, 0) {

      // Check the size of the input
      DIALS_ASSERT(coord.size() > 0);
      DIALS_ASSERT(coord.size() == value.size());

      // Calculate the centroid and total counts
      for (std::size_t i = 0; i < coord.size(); ++i) {
        int c = value[i];
        position_ += (coord[i] * (double)c);
        counts_ += c;
      }

      DIALS_ASSERT(counts_ > 0);
      position_ = position_ / counts_;

      // Calculate the variance on the centroid.
      for (std::size_t i = 0; i < coord.size(); ++i) {
        int c = value[i];
        variance_ += (sqr(coord[i] - position_[0]) * (double)c);
      }
      variance_ = variance_ / counts_;
    }

    /** Get the total counts */
    int counts() {
      return counts_;
    }

    /** Get the centroid position */
    vec2<double> position() {
      return position_;
    }

    /** Get the centroid variance */
    vec2<double> variance() {
      return variance_;
    }

    /** Get the centroid variance per count */
    vec2<double> variance_per_count() {
      DIALS_ASSERT(counts_ > 0);
      return variance_ / counts_ + 1.0;
    }

  private:

    int counts_;
    vec2<double> position_;
    vec2<double> variance_;
  };

  /**
   * Class to calculate the 3D centroid of a list of coordinates
   */
  class CentroidList3d {
  public:

    /**
     * Calculate the centroid.
     * @param coord The list of coordinates.
     * @param value The list of values
     */
    CentroidList3d(const flex_int &value, const flex_vec3_double &coord)
      : counts_(0),
        position_(0, 0, 0),
        variance_(0, 0, 0) {

      // Check the size of the input
      DIALS_ASSERT(coord.size() > 0);
      DIALS_ASSERT(coord.size() == value.size());

      // Calculate the centroid and total counts
      for (std::size_t i = 0; i < coord.size(); ++i) {
        int c = value[i];
        position_ += (coord[i] * (double)c);
        counts_ += c;
      }

      DIALS_ASSERT(counts_ > 0);
      position_ = position_ / counts_;

      // Calculate the variance on the centroid.
      for (std::size_t i = 0; i < coord.size(); ++i) {
        int c = value[i];
        variance_ += (sqr(coord[i] - position_[0]) * (double)c);
      }
      variance_ = variance_ / counts_;
    }

    /** Get the total counts */
    int counts() {
      return counts_;
    }

    /** Get the centroid position */
    vec3<double> position() {
      return position_;
    }

    /** Get the centroid variance */
    vec3<double> variance() {
      return variance_;
    }

    /** Get the centroid variance per count */
    vec3<double> variance_per_count() {
      DIALS_ASSERT(counts_ > 0);
      return variance_ / counts_ + 1.0;
    }

  private:

    int counts_;
    vec3<double> position_;
    vec3<double> variance_;
  };
}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H */
