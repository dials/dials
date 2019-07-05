/*
 * sort_index.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_SORT_INDEX_H
#define DIALS_ARRAY_FAMILY_SORT_INDEX_H

#include <algorithm>

namespace dials { namespace af {

  /**
   * Functor to compare in sort_index.
   */
  template <class RandomAccessIterator>
  struct index_less {
    index_less(const RandomAccessIterator& v) : v_(v) {}

    template <class IndexType>
    bool operator()(const IndexType& x, const IndexType& y) const {
      return v_[x] < v_[y];
    }
    const RandomAccessIterator& v_;
  };

  /**
   * Given a vector return a sorted list of indices.
   * @param v The list of values
   * @returns A sorted list of indices
   */
  template <typename IndexIterator, typename RandomAccessIterator>
  void sort_index(IndexIterator begin, IndexIterator end, RandomAccessIterator v) {
    std::sort(begin, end, index_less<RandomAccessIterator>(v));
  }

}}  // namespace dials::af

#endif /* DIALS_ARRAY_FAMILY_SORT_INDEX_H */
