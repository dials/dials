/*
 * binner.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_BINNER_H
#define DIALS_ARRAY_FAMILY_BINNER_H

#include <map>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace af {

  /**
   * A class to compute the count, sum and mean of values in bins
   */
  class BinIndexer {
  public:
    /**
     * Initialize with the number of bins and an index
     */
    BinIndexer(std::size_t nbins, af::shared<std::size_t> index)
        : nbins_(nbins), index_(index) {
      for (std::size_t i = 0; i < index_.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
      }
    }

    /**
     * Get the indices for items in a bin
     * @param index The bin number
     * @return The indices
     */
    af::shared<std::size_t> indices(std::size_t index) {
      DIALS_ASSERT(index < nbins_);
      af::shared<std::size_t> result;
      for (std::size_t i = 0; i < index_.size(); ++i) {
        if (index_[i] == index) {
          result.push_back(i);
        }
      }
      return result;
    }

    /**
     * @returns A count of the values in each bin
     */
    af::shared<std::size_t> count() const {
      af::shared<std::size_t> result(nbins_, 0);
      for (std::size_t i = 0; i < index_.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
        result[index_[i]]++;
      }
      return result;
    }

    /**
     * @param y The quantity
     * @returns The sum of y in each bin
     */
    af::shared<double> sum(const af::const_ref<double> &y) const {
      DIALS_ASSERT(y.size() == index_.size());
      af::shared<double> result(nbins_, 0);
      for (std::size_t i = 0; i < y.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
        result[index_[i]] += y[i];
      }
      return result;
    }

    /**
     * @param y The quantity
     * @returns The sum of y in each bin
     */
    af::shared<int> sum(const af::const_ref<int> &y) const {
      DIALS_ASSERT(y.size() == index_.size());
      af::shared<int> result(nbins_, 0);
      for (std::size_t i = 0; i < y.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
        result[index_[i]] += y[i];
      }
      return result;
    }

    /**
     * @param y The quantity
     * @returns The sum of y in each bin
     */
    af::shared<int> sum(const af::const_ref<bool> &y) const {
      DIALS_ASSERT(y.size() == index_.size());
      af::shared<int> result(nbins_, 0);
      for (std::size_t i = 0; i < y.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
        result[index_[i]] += (int)y[i];
      }
      return result;
    }

    /**
     * @param y The quantity
     * @returns The mean of y in each bin
     */
    af::shared<double> mean(const af::const_ref<double> &y) const {
      DIALS_ASSERT(y.size() == index_.size());
      af::shared<std::size_t> num = count();
      af::shared<double> result(nbins_, 0);
      for (std::size_t i = 0; i < y.size(); ++i) {
        DIALS_ASSERT(index_[i] < nbins_);
        result[index_[i]] += y[i];
      }
      for (std::size_t i = 0; i < result.size(); ++i) {
        if (num[i] > 0) {
          result[i] /= num[i];
        }
      }
      return result;
    }

  private:
    std::size_t nbins_;
    af::shared<std::size_t> index_;
  };

  /**
   * A class to help with binned data
   */
  class Binner {
  public:
    typedef std::map<double, std::size_t> map_type;

    Binner(const af::const_ref<double> &bins) {
      DIALS_ASSERT(bins.size() > 1);
      bins_.insert(std::pair<double, std::size_t>(bins[0], 0));
      for (std::size_t i = 1; i < bins.size(); ++i) {
        DIALS_ASSERT(bins[i] > bins[i - 1]);
        bins_.insert(bins_.end(), std::pair<double, std::size_t>(bins[i], i));
      }
    }

    /**
     * @returns The bins
     */
    af::shared<double> bins() const {
      af::shared<double> result(bins_.size());
      std::size_t i = 0;
      for (map_type::const_iterator it = bins_.begin(); it != bins_.end(); ++it) {
        DIALS_ASSERT(it->second == i++);
        result[it->second] = it->first;
      }
      return result;
    }

    /**
     * @param x The x value
     * @returns an indexer
     */
    BinIndexer indexer(const af::const_ref<double> &x) {
      // Find the indices of elements
      af::shared<std::size_t> index(x.size());
      for (std::size_t i = 0; i < x.size(); ++i) {
        map_type::const_iterator it = bins_.upper_bound(x[i]);
        if (it != bins_.begin()) {
          --it;
        }
        index[i] = it->second;
      }

      // Return the indexer
      return BinIndexer(bins_.size(), index);
    }

    /**
     * @returns The number of bins
     */
    std::size_t size() const {
      return bins_.size();
    }

  private:
    map_type bins_;
  };

}}  // namespace dials::af

#endif  // DIALS_ARRAY_FAMILY_BINNER_H
