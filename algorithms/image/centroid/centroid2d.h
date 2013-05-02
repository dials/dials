/*
 * centroid2d.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID2D_H
#define DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID2D_H

#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::af::flex_int;

  template <typename T>
  T sqr(const T &a) {
    return a * a;
  }

  /**
   * Class to calculate the 2D centroid of an image.
   */
  class Centroid2d {
  public:

    /**
     * Calculate the centroid.
     * @param image The image to calculate the centroid from.
     */
    Centroid2d(const flex_int &image)
      : counts_(0),
        position_(0, 0),
        variance_(0, 0) {

      // Check the image size
      DIALS_ASSERT(image.accessor().all().size() == 2);
      std::size_t ysize = image.accessor().all()[0];
      std::size_t xsize = image.accessor().all()[1];
      DIALS_ASSERT(xsize > 0 && ysize > 0);

      // Calculate the centroid and total counts
      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {
          int c = image(j, i);
          position_[0] += c * (i + 0.5);
          position_[1] += c * (j + 0.5);
          counts_ += c;
        }
      }

      DIALS_ASSERT(counts_ > 0);
      position_ = position_ / counts_;

      // Calculate the variance on the centroid.
      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {
          int c = image(j, i);
          variance_[0] += c * sqr(i + 0.5 - position_[0]);
          variance_[1] += c * sqr(j + 0.5 - position_[1]);
        }
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
   * Class to calculate the masked 2D centroid of an image.
   */
  class MaskedCentroid2d {
  public:

    /**
     * Calculate the centroid.
     * @param image The image to calculate the centroid from.
     * @param mask The image mask
     */
    MaskedCentroid2d(const flex_int &image, const flex_int &mask)
      : counts_(0),
        position_(0, 0),
        variance_(0, 0) {

      // Check the image size
      DIALS_ASSERT(image.accessor().all().size() == 2);
      DIALS_ASSERT(image.accessor().all().all_eq(mask.accessor().all()));
      std::size_t ysize = image.accessor().all()[0];
      std::size_t xsize = image.accessor().all()[1];
      DIALS_ASSERT(xsize > 0 && ysize > 0);

      // Calculate the centroid and total counts
      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {
          if (mask(j, i)) {
            int c = image(j, i);
            position_[0] += c * (i + 0.5);
            position_[1] += c * (j + 0.5);
            counts_ += c;
          }
        }
      }

      DIALS_ASSERT(counts_ > 0);
      position_ = position_ / counts_;

      // Calculate the variance on the centroid.
      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {
          if (mask(j, i)) {
            int c = image(j, i);
            variance_[0] += c * sqr(i + 0.5 - position_[0]);
            variance_[1] += c * sqr(j + 0.5 - position_[1]);
          }
        }
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
}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID2D_H */
