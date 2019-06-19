/*
 * correlation.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_STATISTICS_CORRELATION_H
#define DIALS_ALGORITHMS_STATISTICS_CORRELATION_H

#include <algorithm>
#include <dials/error.h>

namespace dials { namespace algorithms {

  namespace detail {

    // Struct to help sort
    template <typename T>
    struct sort_by_index {
      const af::const_ref<T> &d;
      sort_by_index(const af::const_ref<T> &d_) : d(d_) {}
      bool operator()(std::size_t a, std::size_t b) const {
        return d[a] < d[b];
      }
    };

  }  // namespace detail

  /**
   * A function to compute the rank of an array
   * @param data The data to rank
   * @return The rank
   */
  template <typename T>
  af::shared<T> rank(const af::const_ref<T> data) {
    // Construct the indices
    std::vector<std::size_t> index(data.size());
    for (std::size_t i = 0; i < index.size(); ++i) {
      index[i] = i;
    }

    // Sort by index
    std::sort(index.begin(), index.end(), detail::sort_by_index<T>(data));

    // Loop through
    af::shared<T> result(index.size());
    for (std::size_t i = 0; i < index.size();) {
      std::size_t j = i + 1;
      T value = (T)(i + 1.0);
      for (; j < index.size(); ++j) {
        DIALS_ASSERT(data[index[j]] >= data[index[i]]);
        if (data[index[j]] > data[index[i]]) {
          break;
        }
        value += (T)(j + 1.0);
      }
      value /= (T)(j - i);
      for (; i < j; ++i) {
        result[index[i]] = value;
      }
    }

    // Return the rank
    return result;
  }

  /**
   * Compute the rank correlation between two arrays
   * @param a an array
   * @param b an array
   * @return The rank correlation coefficient
   */
  template <typename T>
  T spearman_correlation_coefficient(const af::const_ref<T> &a,
                                     const af::const_ref<T> &b) {
    DIALS_ASSERT(a.size() == b.size());

    // Rank the two datasets
    af::shared<T> ra = rank(a);
    af::shared<T> rb = rank(b);

    // The numerator
    T num = 0.0;
    for (std::size_t i = 0; i < ra.size(); ++i) {
      num += (ra[i] - rb[i]) * (ra[i] - rb[i]);
    }
    num *= 6.0;

    // The denominator
    T n = (T)ra.size();
    T den = n * (n * n - 1.0);
    DIALS_ASSERT(den > 0);

    // Return the correlation
    return 1.0 - num / den;
  }

  /**
   * Compute the correlation coefficient between two arrays
   * @param a an array
   * @param b an array
   * @return The corrleation coefficient
   */
  template <typename T>
  T pearson_correlation_coefficient(const af::const_ref<T> &x,
                                    const af::const_ref<T> &y) {
    DIALS_ASSERT(x.size() == y.size());
    DIALS_ASSERT(x.size() > 0);
    T mx = 0.0;
    T my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      mx += x[i];
      my += y[i];
    }
    mx /= x.size();
    my /= y.size();
    T sdx2 = 0.0;
    T sdy2 = 0.0;
    T sdxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      T dx = x[i] - mx;
      T dy = y[i] - my;
      sdx2 += dx * dx;
      sdy2 += dy * dy;
      sdxy += dx * dy;
    }
    DIALS_ASSERT(sdx2 > 0 && sdy2 > 0);
    return sdxy / (std::sqrt(sdx2) * std::sqrt(sdy2));
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_CORRELATION_H
