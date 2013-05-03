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

  using scitbx::af::flex;
  using scitbx::af::flex_int;

  template <typename T>
  T sqr(const T &a) {
    return a * a;
  }

  /**
   * Class to calculate the centroid of a list of coordinates
   */
  template <typename CoordType>
  class CentroidList {
  public:

    // Useful typedefs
    typedef CoordType coord_type;
    typedef typename flex<coord_type>::type flex_type;

    /**
     * Calculate the centroid.
     * @param coord The list of coordinates.
     * @param value The list of values
     */
    CentroidList(const flex_int &value, const flex_type &coord)
      : counts_(0),
        position_(0.0),
        squared_width_(0.0) {

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
        squared_width_ += (sqr(coord[i] - position_[0]) * (double)c);
      }
      squared_width_ = squared_width_ / counts_;
    }

    /** Get the total counts */
    int counts() {
      return counts_;
    }

    /** Get the centroid position */
    coord_type position() {
      return position_;
    }

    /** Get the centroid variance */
    coord_type squared_width() {
      return squared_width_;
    }

    /** Get the centroid variance per count */
    coord_type variance() {
      DIALS_ASSERT(counts_ > 0);
      return squared_width_ / counts_ + 1.0;
    }

  private:

    int counts_;
    coord_type position_;
    coord_type squared_width_;
  };
}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H */
