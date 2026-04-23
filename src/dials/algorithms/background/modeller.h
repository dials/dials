/*
 * modeller.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_MODELLER_H
#define DIALS_ALGORITHMS_BACKGROUND_MODELLER_H

#include <dials/model/data/image_volume.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace dials { namespace algorithms {

  using dials::model::Foreground;
  using dials::model::ImageVolume;
  using dials::model::MultiPanelImageVolume;
  using dials::model::Valid;

  template <typename ElementType>
  class af_shared_with_getitem : public af::shared<ElementType> {
  public:
    ElementType getitem(const int i) {
      return (*this)[i];
    }
  };

  /**
   * Class to compute background statistics
   */
  class BackgroundStatistics {
  public:
    /**
     * Initialize from an image volume
     * @param volume The image volume
     */
    BackgroundStatistics(const ImageVolume<>& volume)
        : accessor_(volume.accessor()[1], volume.accessor()[2]),
          sum_(accessor_, 0.0),
          sum_sq_(accessor_, 0.0),
          num_(accessor_, 0),
          min_(accessor_, -1),
          max_(accessor_, -1) {
      DIALS_ASSERT(volume.is_consistent());
      typedef ImageVolume<>::float_type FloatType;
      af::const_ref<FloatType, af::c_grid<3> > data = volume.data().const_ref();
      af::const_ref<int, af::c_grid<3> > mask = volume.mask().const_ref();
      for (std::size_t j = 0; j < accessor_[0]; ++j) {
        for (std::size_t i = 0; i < accessor_[1]; ++i) {
          for (std::size_t k = 0; k < data.accessor()[0]; ++k) {
            double d = data(k, j, i);
            int m = mask(k, j, i);
            if ((m & Valid) && !(m & Foreground)) {
              sum_(j, i) += d;
              sum_sq_(j, i) += d * d;
              num_(j, i) += 1;
              if (min_(j, i) == -1 || min_(j, i) > d) min_(j, i) = d;
              if (max_(j, i) == -1 || max_(j, i) < d) max_(j, i) = d;
            }
          }
        }
      }
    }

    p::object get_all() const {
      long acc0 = accessor_[0];
      long acc1 = accessor_[1];

      p::tuple shape = p::make_tuple(acc0, acc1);
      p::tuple stride_double = p::make_tuple(acc1 * sizeof(double), sizeof(double));
      p::tuple stride_int = p::make_tuple(acc1 * sizeof(int), sizeof(int));
      np::ndarray sum = np::from_data(sum_.begin(),
                                      np::dtype::get_builtin<double>(),
                                      shape,
                                      stride_double,
                                      p::object());
      np::ndarray sum_sq = np::from_data(sum_sq_.begin(),
                                         np::dtype::get_builtin<double>(),
                                         shape,
                                         stride_double,
                                         p::object());
      np::ndarray num = np::from_data(
        num_.begin(), np::dtype::get_builtin<int>(), shape, stride_int, p::object());
      np::ndarray min = np::from_data(min_.begin(),
                                      np::dtype::get_builtin<double>(),
                                      shape,
                                      stride_double,
                                      p::object());
      np::ndarray max = np::from_data(max_.begin(),
                                      np::dtype::get_builtin<double>(),
                                      shape,
                                      stride_double,
                                      p::object());

      // np::ndarray sum = np::empty(2, shape, np::dtype::get_builtin<double>());
      // np::ndarray sum_sq = np::empty(2, shape, np::dtype::get_builtin<double>());
      // np::ndarray num = np::empty(2, shape, np::dtype::get_builtin<int>());
      // np::ndarray min = np::empty(2, shape, np::dtype::get_builtin<double>());
      // np::ndarray max = np::empty(2, shape, np::dtype::get_builtin<double>());

      // double* s = reinterpret_cast<double*>(sum.get_data());
      // double* ss = reinterpret_cast<double*>(sum_sq.get_data());
      // int* n = reinterpret_cast<int*>(num.get_data());
      // double* mi = reinterpret_cast<double*>(min.get_data());
      // double* ma = reinterpret_cast<double*>(max.get_data());

      // for (int i = 0; i < acc0; i++) {
      //   for (int j = 0; j < acc1; j++) {
      //     s[i * acc1 + j] = sum_(i, j);
      //     ss[i * acc1 + j] = sum_sq_(i, j);
      //     n[i * acc1 + j] = num_(i, j);
      //     mi[i * acc1 + j] = min_(i, j);
      //     ma[i * acc1 + j] = max_(i, j);
      //   }
      // }

      return p::make_tuple(
        acc0, acc1, sum.copy(), sum_sq.copy(), num.copy(), min.copy(), max.copy());
    }

    void set_all(p::object state) {
      const long acc0 = p::extract<long>(state[0]);
      const long acc1 = p::extract<long>(state[1]);
      const np::ndarray sum = p::extract<np::ndarray>(state[2]);
      const np::ndarray sum_sq = p::extract<np::ndarray>(state[3]);
      const np::ndarray num = p::extract<np::ndarray>(state[4]);
      const np::ndarray min = p::extract<np::ndarray>(state[5]);
      const np::ndarray max = p::extract<np::ndarray>(state[6]);

      double* s = reinterpret_cast<double*>(sum.get_data());
      double* ss = reinterpret_cast<double*>(sum_sq.get_data());
      int* n = reinterpret_cast<int*>(num.get_data());
      double* mi = reinterpret_cast<double*>(min.get_data());
      double* ma = reinterpret_cast<double*>(max.get_data());

      accessor_ = af::c_grid<2>(acc0, acc1);
      sum_ = af::versa<double, af::c_grid<2> >(accessor_);
      sum_sq_ = af::versa<double, af::c_grid<2> >(accessor_);
      num_ = af::versa<int, af::c_grid<2> >(accessor_);
      min_ = af::versa<double, af::c_grid<2> >(accessor_);
      max_ = af::versa<double, af::c_grid<2> >(accessor_);

      size_t n_elems = static_cast<size_t>(acc0) * static_cast<size_t>(acc1);

      std::memcpy(sum_.begin(), s, n_elems * sizeof(double));
      std::memcpy(sum_sq_.begin(), ss, n_elems * sizeof(double));
      std::memcpy(num_.begin(), n, n_elems * sizeof(int));
      std::memcpy(min_.begin(), mi, n_elems * sizeof(double));
      std::memcpy(max_.begin(), ma, n_elems * sizeof(double));
    }

    /**
     * Add results from another object
     * @param other The other object
     */
    BackgroundStatistics operator+=(const BackgroundStatistics& other) {
      DIALS_ASSERT(accessor_.all_eq(other.accessor_));
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        sum_[i] += other.sum_[i];
        sum_sq_[i] += other.sum_sq_[i];
        num_[i] += other.num_[i];
        if (min_[i] == -1) {
          min_[i] = other.min_[i];
        } else if (other.min_[i] != -1) {
          min_[i] = std::min(min_[i], other.min_[i]);
        }
        if (max_[i] == -1) {
          max_[i] = other.max_[i];
        } else if (other.max_[i] != -1) {
          max_[i] = std::max(max_[i], other.max_[i]);
        }
      }
      return *this;
    }

    /**
     * @returns The image sum at each pixel
     */
    af::versa<double, af::c_grid<2> > sum() const {
      return sum_;
    }

    /**
     * @returns The image sum_sq at each pixel
     */
    af::versa<double, af::c_grid<2> > sum_sq() const {
      return sum_sq_;
    }

    /**
     * @returns The number of images contributing for each pixel
     */
    af::versa<int, af::c_grid<2> > num() const {
      return num_;
    }

    /**
     * @returns The minimum image
     */
    af::versa<double, af::c_grid<2> > min() const {
      return min_;
    }

    /**
     * @returns The maximum image
     */
    af::versa<double, af::c_grid<2> > max() const {
      return max_;
    }

    /**
     * @returns The mean at each pixel
     */
    af::versa<double, af::c_grid<2> > mean(std::size_t min_images) const {
      DIALS_ASSERT(min_images > 0);
      af::versa<double, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (num_[i] >= min_images) {
          result[i] = sum_[i] / num_[i];
        }
      }
      return result;
    }

    /**
     * @returns The variance at each pixel
     */
    af::versa<double, af::c_grid<2> > variance(std::size_t min_images) const {
      DIALS_ASSERT(min_images > 0);
      af::versa<double, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (num_[i] >= min_images) {
          result[i] = (sum_sq_[i] - sum_[i] * sum_[i] / num_[i]) / num_[i];

          // Sometimes this becomes < 0 due to numerical instability
          if (result[i] < 0) {
            result[i] = 0;
          }

          DIALS_ASSERT(result[i] >= 0);
        }
      }
      return result;
    }

    /**
     * @returns The dispersion at each pixel
     */
    af::versa<double, af::c_grid<2> > dispersion(std::size_t min_images) const {
      DIALS_ASSERT(min_images > 0);
      af::versa<double, af::c_grid<2> > m = mean(min_images);
      af::versa<double, af::c_grid<2> > v = variance(min_images);
      af::versa<double, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (m[i] > 0) {
          result[i] = v[i] / m[i];
        }
      }
      return result;
    }

    /**
     * @returns The mask at each pixel
     */
    af::versa<bool, af::c_grid<2> > mask(std::size_t min_images) const {
      DIALS_ASSERT(min_images > 0);
      af::versa<bool, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = num_[i] >= min_images;
      }
      return result;
    }

  private:
    af::c_grid<2> accessor_;
    af::versa<double, af::c_grid<2> > sum_;
    af::versa<double, af::c_grid<2> > sum_sq_;
    af::versa<int, af::c_grid<2> > num_;
    af::versa<double, af::c_grid<2> > min_;
    af::versa<double, af::c_grid<2> > max_;
  };

  /**
   * A class to do multi panel background statistics
   */
  class MultiPanelBackgroundStatistics {
  public:
    /**
     * Initialize with multipanel image volume
     * @param volume The multi panel image volume
     */
    MultiPanelBackgroundStatistics(const MultiPanelImageVolume<>& volume) {
      for (std::size_t i = 0; i < volume.size(); ++i) {
        statistics_.push_back(BackgroundStatistics(volume.get(i)));
      }
    }

    af_shared_with_getitem<BackgroundStatistics> get_statistics() const {
      return statistics_;
    }

    void set_statistics(af_shared_with_getitem<BackgroundStatistics> statistics) {
      statistics_ = statistics;
    }

    /**
     * @returns the statistics for the given panel
     */
    BackgroundStatistics get(std::size_t index) const {
      DIALS_ASSERT(index < statistics_.size());
      return statistics_[index];
    }

    /**
     * @returns The number of panels
     */
    std::size_t size() const {
      return statistics_.size();
    }

    /**
     * Add results from another object
     * @param other The other object
     */
    MultiPanelBackgroundStatistics operator+=(
      const MultiPanelBackgroundStatistics& other) {
      DIALS_ASSERT(size() == other.size());
      for (std::size_t i = 0; i < size(); ++i) {
        statistics_[i] += other.statistics_[i];
      }
      return *this;
    }

  private:
    af_shared_with_getitem<BackgroundStatistics> statistics_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_BACKGROUND_MODELLER_H
