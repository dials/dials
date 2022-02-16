#ifndef DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H
#define DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H

#include <scitbx/array_family/shared.h>
#include <dials/error.h>


namespace dials { namespace algorithms {

  template<typename T>
  static inline double Lerp(T v0, T v1, T t)
  {
    return (1 - t)*v0 + t*v1;
  }

  class BinnedStatistics
  {
  public:
    BinnedStatistics(scitbx::af::const_ref< double > values,
              scitbx::af::const_ref< std::size_t > bins,
              std::size_t n_bins)
      : values_(values),
        bins_(bins),
        n_bins_(n_bins)
       {
          DIALS_ASSERT(values.size() == bins.size());
          binned_values.reserve(n_bins_);
          is_sorted.reserve(n_bins_);
          for (std::size_t i = 0; i < n_bins; i++)
          {
            binned_values.push_back(scitbx::af::shared<double>());
            is_sorted.push_back(false);
          }
          for (std::size_t i = 0; i < values.size(); i++)
          {
            binned_values[bins_[i]].push_back(values_[i]);
          }
       };

    scitbx::af::shared<double> get_values_in_bin(std::size_t i) const
    {
      DIALS_ASSERT(i < n_bins_);
      return binned_values[i];
    };

    scitbx::af::shared<bool> bin_is_empty() const
    {
      scitbx::af::shared<bool> empty;
      for (std::size_t i = 0; i < n_bins_; i++)
      {
        empty.push_back(binned_values[i].size() == 0);
      }
      return empty;
    };

    scitbx::af::shared<bool> bin_is_sorted() const {
      return is_sorted;
    };

    scitbx::af::shared<double> get_medians()
    {
      scitbx::af::shared<double> median;
      for (std::size_t i = 0; i < n_bins_; i++)
      {
        if (not is_sorted[i])
        {
          std::sort(binned_values[i].begin(), binned_values[i].end());
          is_sorted[i] = true;
        }
        std::size_t n = binned_values[i].size();

        if (n == 0) median.push_back(0);
        else if (n == 1) median.push_back(binned_values[i][0]);
        else if (n == 2) median.push_back((binned_values[i][0] + binned_values[i][1])/2);
        else
        {
          double * const k = binned_values[i].begin() + n/2;
          if (n % 2) median.push_back(*k);
          else
          {
            double * const k1 = binned_values[i].begin() + n/2 - 1;
            median.push_back((*k + *k1)/2);
          }
        }
      }
      return median;
    };

    scitbx::af::shared<double> get_iqrs()
    {
      // https://stackoverflow.com/questions/11964552/finding-quartiles
      scitbx::af::shared<double> iqr;
      for (std::size_t i = 0; i < n_bins_; i++)
      {
        if (not is_sorted[i])
        {
          std::sort(binned_values[i].begin(), binned_values[i].end());
          is_sorted[i] = true;
        }
        std::size_t n = binned_values[i].size();

        if (n == 0) iqr.push_back(0);
        else if (n == 1) iqr.push_back(0);
        else
        {
          // q1
          double poi = Lerp<double>(-0.5, n -0.5, 0.25);
          size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
          size_t right = std::min(int64_t(std::ceil(poi)), int64_t(n - 1));
          double dat_left = binned_values[i].at(left);
          double dat_right = binned_values[i].at(right);
          double q1 = Lerp<double>(dat_left, dat_right, poi - left);

          // q3
          poi = Lerp<double>(-0.5, n -0.5, 0.75);
          left = std::max(int64_t(std::floor(poi)), int64_t(0));
          right = std::min(int64_t(std::ceil(poi)), int64_t(n - 1));
          dat_left = binned_values[i].at(left);
          dat_right = binned_values[i].at(right);
          double q3 = Lerp<double>(dat_left, dat_right, poi - left);

          iqr.push_back(q3 - q1);
        }
      }
      return iqr;
    };

  private:
    scitbx::af::const_ref< double > values_;
    scitbx::af::const_ref< std::size_t > bins_;
    std::size_t n_bins_;

    std::vector<scitbx::af::shared<double> > binned_values;
    scitbx::af::shared<bool> is_sorted;
  };


}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H
