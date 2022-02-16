#ifndef DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H
#define DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H

#include <scitbx/array_family/shared.h>
//#include <scitbx/array_family/ref.h>
//#include <scitbx/array_family/versa.h>
#include <dials/error.h>


namespace dials { namespace algorithms {

  class BinnedStatistics {
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
          for (std::size_t i = 0; i < n_bins; i++){
            binned_values.push_back(scitbx::af::shared<double>());
          }
          for (std::size_t i = 0; i < values.size(); i++){
            binned_values[bins_[i]].push_back(values_[i]);
          }
       }

    scitbx::af::shared<double> get_values_in_bin(std::size_t i) const {
      DIALS_ASSERT(i < n_bins_);
      return binned_values[i];
    };

    scitbx::af::shared<bool> bin_is_empty() const {
      scitbx::af::shared<bool> empty;
      for (std::size_t i = 0; i < n_bins_; i++){
        empty.push_back(binned_values[i].size() == 0);
      }
      return empty;
    }

  private:
    scitbx::af::const_ref< double > values_;
    scitbx::af::const_ref< std::size_t > bins_;
    std::size_t n_bins_;

    std::vector<scitbx::af::shared<double> > binned_values;
  };


}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_BINNED_STATISTICS_H
