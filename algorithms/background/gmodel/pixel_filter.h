/*
 * pixel_filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_GMODEL_PIXEL_FILTER_H
#define DIALS_ALGORITHMS_BACKGROUND_GMODEL_PIXEL_FILTER_H

#include <vector>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  /**
   * A class to hold the value from filtering
   */
  class PixelFilterResult {
  public:
    /**
     * Initialise the class
     * @param data The data array
     * @param mask The mask array
     */
    PixelFilterResult(af::versa<double, af::c_grid<2> > data,
                      af::versa<bool, af::c_grid<2> > mask)
        : data_(data), mask_(mask) {
      DIALS_ASSERT(data_.accessor().all_eq(mask_.accessor()));
    }

    /**
     * @return The data array
     */
    af::versa<double, af::c_grid<2> > data() const {
      return data_;
    }

    /**
     * @return The mask array
     */
    af::versa<bool, af::c_grid<2> > mask() const {
      return mask_;
    }

  private:
    af::versa<double, af::c_grid<2> > data_;
    af::versa<bool, af::c_grid<2> > mask_;
  };

  /**
   * A class to filter the image pixels
   */
  class PixelFilter {
  public:
    /**
     * Initialise the pixel filter
     * @param height The height of the image
     * @param width The width of the image
     */
    PixelFilter(std::size_t height, std::size_t width)
        : max_count_(0),
          height_(height),
          width_(width),
          sum1_(height * width),
          sum2_(height * width),
          count_(height * width) {}

    /**
     * Add another image to be processed
     * @param data The image data
     * @param mask The image mask
     */
    template <typename T>
    void add(const af::const_ref<T, af::c_grid<2> > &data,
             const af::const_ref<bool, af::c_grid<2> > &mask) {
      // Check stuff is ok
      DIALS_ASSERT(data.accessor()[0] == height_);
      DIALS_ASSERT(data.accessor()[1] == width_);
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

      // Update the sums
      for (std::size_t i = 0; i < data.size(); ++i) {
        if (mask[i]) {
          count_[i] += 1;
          sum1_[i] += (double)data[i];
          sum2_[i] += (double)data[i] * (double)data[i];
        }
      }

      // Increment max count
      max_count_++;
    }

    /**
     * Compute the filtered image
     * @param min_count The minimum number of elements per pixel
     * @param nsigma The number of standard deviations to filter by
     */
    PixelFilterResult compute(std::size_t min_count, double nsigma) const {
      // Check input
      DIALS_ASSERT(nsigma > 0);
      DIALS_ASSERT(max_count_ >= 2);
      if (min_count < 2 || min_count > max_count_) {
        min_count = max_count_;
      }

      // Initialise the result array
      af::versa<double, af::c_grid<2> > data(af::c_grid<2>(height_, width_), 0);
      af::versa<bool, af::c_grid<2> > mask(af::c_grid<2>(height_, width_), false);

      // Pixels whose variance is > expected variance are masked out using the
      // index of dispersion. Otherwise, the result is the mean value
      for (std::size_t i = 0; i < count_.size(); ++i) {
        if (count_[i] >= min_count) {
          double s1 = sum1_[i];
          double s2 = sum2_[i];
          double n = count_[i];
          double mean = s1 / n;
          double var = (s2 - s1 * s1 / n) / (n - 1);
          if (var <= mean * (1.0 + nsigma * sqrt(2.0 / (n - 1)))) {
            data[i] = mean;
            mask[i] = true;
          }
        }
      }

      // Return the pixel filter value
      return PixelFilterResult(data, mask);
    }

    /**
     * @return The number of images
     */
    std::size_t num_images() const {
      return max_count_;
    }

  private:
    std::size_t max_count_;
    std::size_t height_;
    std::size_t width_;
    std::vector<double> sum1_;
    std::vector<double> sum2_;
    std::vector<std::size_t> count_;
  };

}}}  // namespace dials::algorithms::background

#endif  // DIALS_ALGORITHMS_BACKGROUND_GMODEL_PIXEL_FILTER_H
