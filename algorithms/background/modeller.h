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

namespace dials { namespace algorithms {

  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::model::ImageVolume;
  using dials::model::MultiPanelImageVolume;

  /**
   * Class to compute background statistics
   */
  class BackgroundStatistics {
  public:

    /**
     * Initialize from an image volume
     * @param volume The image volume
     */
    BackgroundStatistics(const ImageVolume<> &volume)
      : accessor_(
          volume.accessor()[1],
          volume.accessor()[2]),
        sum_(accessor_, 0.0),
        sum_sq_(accessor_, 0.0),
        num_(accessor_, 0) {
      DIALS_ASSERT(volume.is_consistent());
      typedef ImageVolume<>::float_type FloatType;
      af::const_ref< FloatType, af::c_grid<3> > data = volume.data().const_ref();
      af::const_ref< int, af::c_grid<3> > mask = volume.mask().const_ref();
      for (std::size_t j = 0; j < accessor_[0]; ++j) {
        for (std::size_t i = 0; i < accessor_[1]; ++i) {
          for (std::size_t k = 0; k < data.accessor()[0]; ++k) {
            double d = data(k,j,i);
            int m = mask(k,j,i);
            if ((m & Valid) && !(m & Foreground)) {
              sum_(j,i) += d;
              sum_sq_(j,i) += d * d;
              num_(j,i) += 1;
            }
          }
        }
      }
    }

    /**
     * Add results from another object
     * @param other The other object
     */
    BackgroundStatistics operator+=(const BackgroundStatistics &other) {
      DIALS_ASSERT(accessor_.all_eq(other.accessor_));
      for (std::size_t i = 0; i < sum_.size(); ++i) {
        sum_[i] += other.sum_[i];
        sum_sq_[i] += other.sum_sq_[i];
        num_[i] += other.num_[i];
      }
      return *this;
    }

    /**
     * @returns The image sum at each pixel
     */
    af::versa< double, af::c_grid<2> > sum() const {
      return sum_;
    }

    /**
     * @returns The image sum_sq at each pixel
     */
    af::versa< double, af::c_grid<2> > sum_sq() const {
      return sum_sq_;
    }

    /**
     * @returns The number of images contributing for each pixel
     */
    af::versa< int, af::c_grid<2> > num() const {
      return num_;
    }

    /**
     * @returns The mean at each pixel
     */
    af::versa< double, af::c_grid<2> > mean() const {
      af::versa <double, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (num_[i] > 0) {
          result[i] = sum_[i] / num_[i];
        }
      }
      return result;
    }

    /**
     * @returns The variance at each pixel
     */
    af::versa< double, af::c_grid<2> > variance() const {
      af::versa <double, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (num_[i] > 0) {
          result[i] = (sum_sq_[i] - sum_[i]*sum_[i] / num_[i]) / num_[i];
          DIALS_ASSERT(result[i] >= 0);
        }
      }
      return result;
    }

    /**
     * @returns The dispersion at each pixel
     */
    af::versa< double, af::c_grid<2> > dispersion() const {
      af::versa<double, af::c_grid<2> > m = mean();
      af::versa<double, af::c_grid<2> > v = variance();
      af::versa <double, af::c_grid<2> > result(accessor_);
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
    af::versa < bool, af::c_grid<2> > mask() const {
      af::versa <bool, af::c_grid<2> > result(accessor_);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = num_[i] > 0;
      }
      return result;
    }

  private:

    af::c_grid<2> accessor_;
    af::versa< double, af::c_grid<2> > sum_;
    af::versa< double, af::c_grid<2> > sum_sq_;
    af::versa< int, af::c_grid<2> > num_;
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
    MultiPanelBackgroundStatistics(const MultiPanelImageVolume<> &volume) {
      for (std::size_t i = 0; i < volume.size(); ++i) {
        statistics_.push_back(BackgroundStatistics(volume.get(i)));
      }
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
    MultiPanelBackgroundStatistics operator+=(const MultiPanelBackgroundStatistics &other) {
      DIALS_ASSERT(size() == other.size());
      for (std::size_t i = 0; i < size(); ++i) {
        statistics_[i] += other.statistics_[i];
      }
      return *this;
    }

  private:

    af::shared<BackgroundStatistics> statistics_;
  };


}}

#endif // DIALS_ALGORITHMS_BACKGROUND_MODELLER_H
