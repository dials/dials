/*
 * detect_collisions.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPATIAL_INDEXING_DETECT_COLLISIONS_H
#define DIALS_ALGORITHMS_SPATIAL_INDEXING_DETECT_COLLISIONS_H

#include <cstdlib>
#include <vector>
#include <algorithm>

namespace dials { namespace algorithms {

  /** Dumb limits for the algorithm to stop recursion */
  enum {
    QC_MAXDEPTH = 20,           ///< Maximum recursion depth
    BF_THRESHOLD = 10           ///< Threshold for brute force search
  };

  /**
   * A simple bounding box structure that has min/max for each
   * dimension 0 <= D < DIM
   *
   * @tparam DIM the number of dimensions
   */
  template <int DIM>
  struct bounding_box {
    double min[DIM];
    double max[DIM];
    bounding_box() {}
    bounding_box(const bounding_box &box) {
      for (std::size_t i = 0; i < DIM; ++i) {
        min[i] = box.min[i];
        max[i] = box.max[i];
      }
    }
  };

  /** Get the minimum bound of the box of type T in dimension DIM */
  template <int DIM, typename T>
  double get_minimum_bound(const T &);

  /** Get the maximum bound of the box of type T in dimension DIM */
  template <int DIM, typename T>
  double get_maximum_bound(const T &);

  /** The init box compile-time loop */
  template <int D, int DIM, typename T>
  struct init_bounding_box_loop {
    static void value(const T &a, bounding_box<DIM> &box) {
      box.min[D-1] = get_minimum_bound<D-1>(a);
      box.max[D-1] = get_maximum_bound<D-1>(a);
      init_bounding_box_loop<D-1, DIM, T>::value(a, box);
    }
  };

  /** Stop the init box compile-time loop */
  template <int DIM, typename T>
  struct init_bounding_box_loop<0, DIM, T> {
    static void value(const T &a, bounding_box<DIM> &box) {}
  };

  /** Initialise the box from the given element */
  template <int DIM, typename T>
  void init_bounding_box(const T &a, bounding_box<DIM> &box) {
    return init_bounding_box_loop<DIM, DIM, T>::value(a, box);
  }

  /** The compile-time loop to set min/max elements */
  template <int D, int DIM, typename T>
  struct set_minmax_bounding_box_loop {
    static void value(const T &a, bounding_box<DIM> &box) {
      if (get_minimum_bound<D-1>(a) < box.min[D-1]) {
        box.min[D-1] = get_minimum_bound<D-1>(a);
      }
      if (get_maximum_bound<D-1>(a) > box.max[D-1]) {
        box.max[D-1] = get_maximum_bound<D-1>(a);
      }
      set_minmax_bounding_box_loop<D-1, DIM, T>::value(a, box);
    }
  };

  /** Stop compile time loop to set min/max elements */
  template <int DIM, typename T>
  struct set_minmax_bounding_box_loop<0, DIM, T> {
    static void value(const T &a, bounding_box<DIM> &box) {}
  };

  /**
   *For the given element, if the lower bounds are less than those
   * of the box, set them, if the upper bounds are greater than then
   * set them.
   */
  template <int DIM, typename T>
  void set_minmax_bounding_box(const T &a, bounding_box<DIM> &box) {
    set_minmax_bounding_box_loop<DIM, DIM, T>::value(a, box);
  }

  /**
   * Get the bounding box from the input data. Basically scan the data for
   * the min/max in each dimension to get the bound that surrounds all the
   * elements. Use some template metaprogramming to cope with different
   * numbers of dimensions.
   */
  template <int DIM, typename Iterator>
  bounding_box<DIM> get_bounding_box(Iterator first, Iterator last) {
    bounding_box<DIM> box;
    init_bounding_box(*first, box);
    for (Iterator it = first+1; it < last; ++it) {
      set_minmax_bounding_box(*it, box);
    }
    return box;
  }

  /**
   * Parition the data based on the lower bound of the data in the given
   * axis. All elements with their lower bound lower than the split to the left,
   * all others to the right.
   */
  template <int DIM, typename Iterator>
  struct partition_by_lower {
    const Iterator &data_;
    double div_;

    partition_by_lower(const Iterator &data, double div)
      : data_(data), div_(div) {}

    bool operator()(const int &a) {
      return get_minimum_bound<DIM>(*(data_ + a)) < div_;
    }
  };

  /**
   * Parition the data based on the upper bound of the data in the given
   * axis. All elements with their upper bound lower than the split to the left,
   * all others to the right.
   */
  template <int DIM, typename Iterator>
  struct partition_by_upper {
    const Iterator &data_;
    double div_;

    partition_by_upper(const Iterator &data, double div)
      : data_(data), div_(div) {}

    bool operator()(const int &a) {
      return get_maximum_bound<DIM>(*(data_ + a)) < div_;
    }
  };

  /**
   * Template stuff to allow creation of struct without explicitly
   * stating type of the iterator we're passing int.
   */
  template <int DIM, typename Iterator>
  partition_by_lower<DIM, Iterator> by_lower(const Iterator &d, double div) {
    return partition_by_lower<DIM, Iterator>(d, div);
  }

  /**
   * Template stuff to allow creation of struct without explicitly
   * stating type of the iterator we're passing int.
   */
  template <int DIM, typename Iterator>
  partition_by_upper<DIM, Iterator> by_upper(const Iterator &d, double div) {
    return partition_by_upper<DIM, Iterator>(d, div);
  }

  /**
   * A simple brute force collision detection algorithm. Each
   * element in the range is checked against every other element. If
   * a collision is encountered then the pair of indices are added
   * to a list. The list assumes there is a method "push_back" where new
   * elements can be added. Additionally it assumes that it defines a
   * dependant type "value_type" (which could be a std::pair<int, int>)
   * whose constructor takes two indices.
   *
   * @tparam Iterator The type of random access iterator
   * @tparam ListType a container that allows items to be appended.
   * @tparam Collides The function that checks if there is a collision.
   *
   * @param first The first iterator in the range
   * @param last The last iterator in the range
   * @param collisions The list of collisions
   * @param collides The collision checking function.
   */
  template <typename Iterator, typename ListType, typename Collides>
  void detect_collisions_brute_force(Iterator first, Iterator last,
      ListType &collisions, Collides collides) {
    // Compare each element against every other. If element a is
    // deemed to collide with element b then add to the list.
    for (Iterator a = first; a < last - 1; ++a) {
      for (Iterator b = a + 1; b < last; ++b) {
        if (collides(*a, *b)) {
          collisions.push_back(typename ListType::value_type(*a, *b));
        }
      }
    }
  }

  /** The compile-time loop for the no_collision test */
  template <int DIM, typename T>
  struct no_collision_loop{
    static bool value(const T &a, const T &b) {
      return get_minimum_bound<DIM-1>(a) > get_maximum_bound<DIM-1>(b)
          || get_minimum_bound<DIM-1>(b) > get_maximum_bound<DIM-1>(a)
          || no_collision_loop<DIM-1, T>::value(a, b);
    }
  };

  /** Stop the compile-time loop for no_collision test */
  template <typename T>
  struct no_collision_loop<0, T>{
    static bool value(const T &a, const T &b) { return false; }
  };

  /**
   * Check if the objects overlap (i.e. in 2D do the following)
   * return !(a.x0 > b.x1 || b.x0 > a.x1 || a.y0 > b.y1 || b.y0 > a.y1);
   *
   * To cope with different numbers of dimensions this is done using a bit
   * of template metaprogramming. The no_collision_loop struct sets up a
   * compile time loop over the number of dimensions.
   */
  template <int DIM, typename T>
  bool collides(const T &a, const T &b) {
    return !no_collision_loop<DIM, T>::value(a, b);
  }

  /** The compile-time loop for already visited test */
  template <int DIM, int BDIM, typename T>
  struct already_visited_loop {
    static bool value(const T &a, const T &b, const bounding_box<BDIM> &box) {
      return (get_minimum_bound<DIM-1>(a) < box.min[DIM-1]
          &&  get_minimum_bound<DIM-1>(b) < box.min[DIM-1])
          || already_visited_loop<DIM-1, BDIM, T>::value(a, b, box);
    }
  };

  /** Stop the compile-time loop for already visited test */
  template <int BDIM, typename T>
  struct already_visited_loop<0, BDIM, T> {
    static bool value(const T &a, const T &b, const bounding_box<BDIM> &box) {
      return false;
    }
  };

  /*
   * If the objects left x coordinate is lower than the split or the
   * objects bottom y coordinate is lower than the split return true.
   * I.e in 2D do the following
   * return (a.x0 < box_.min[0] && b.x0 < box_.min[0])
   *     || (a.y0 < box_.min[1] && b.y0 < box_.min[1]);
   *
   * To cope with different numbers of dimensions this is done using a bit
   * of template metaprogramming. The already_visited_loop struct sets up a
   * compile time loop over the number of dimensions.
   */
  template <int DIM, typename T>
  bool already_visited(const T &a, const T &b, const bounding_box<DIM> &box) {
    return already_visited_loop<DIM, DIM, T>::value(a, b, box);
  }

  /**
   * A struct to do the collision checking as input to the brute-force
   * algorithm. When called, checks if the given objects collide and then
   * performs a further check to see if these objects would have already
   * been added to the collision list. Because of the regular search
   * pattern used by the algorithm, if both colliding objects have their
   * lower bound in a position to the left and down (in 2D) of the current
   * split, then they have already been found and added. Therefore we
   * ignore them.
   */
  template <int DIM, typename Iterator>
  struct check_collision {

    const Iterator d_;
    const bounding_box<DIM> &box_;

    /** Initialise the functor with the iterator and bounding box */
    check_collision(const Iterator &d, const bounding_box<DIM> &box)
      : d_(d), box_(box) {}

    /**
     * Check for a collision. We return true if the objects collide and
     * if they have not already been set.
     */
    template <typename T>
    bool operator()(const T &a, const T &b) {
      // If the objects collide and they have not already been
      // added then return true, otherwise return false.
      return collides<DIM>(*(d_ + a), *(d_ + b))
          && !already_visited<DIM>(*(d_ + a), *(d_ + b), box_);
    }
  };

  /**
   * A helper function to wrap the call to detect_collisions_brute_force with
   * the check_collision struct.
   */
  template <int DIM, typename IndexIterator,
            typename DataIterator, typename ListType>
  void detect_brute_force_w_check(IndexIterator first, IndexIterator last,
      DataIterator data, ListType &collisions, const bounding_box<DIM> &box) {

    // Do a brute force collision test using the collision checker to ensure
    // no pairs are added that have already been visited and added.
    detect_collisions_brute_force(first, last, collisions,
      check_collision<DIM, DataIterator>(data, box));
  }

  /**
   * The main body of the algorithm.
   *
   * This function is templated to partition along each axis (0 <= D < DIM).
   * This allows the same code to be used to instantiate a co-recursive
   * function call. If partition_data is called with D = 0, then the function
   * will then call partition_data with D = 1 which in 2D will call
   * partition_data with D = 0 again etc.
   *
   * The function performs the following tasks. It first checks the exit
   * criteria, namely have we reached the maximum recursion depth or are
   * there fewer elements than the hard coded threshold to do a brute-force
   * search of collisions.
   *
   * If we haven't reached the exit criteria then split the bouding box in
   * half along the current axis. Parition the index array so that all elements
   * with the lower bound of their bounding box less than the split are to the
   * left and all those greater to the right. Then partition along the next
   * axis recursively.
   *
   * Next partition the data such that all elements with their upper bound
   * greater than the split are to the right and all those less than to the
   * less. This will move some elements from the previous partition but is
   * necessary to ensure that those elements that are not wholely within a
   * single partition are properly accounted for. Then partition along the
   * next axis recursively.
   *
   * Once the exit condition has been met, do a brute-force search for
   * collisions amoung the remaining elements.
   *
   * @todo The algorithm uses a dumb max recursion depth to exit, this could
   *    be improved on since it is non-optimal for highly clustered data with
   *    outliers. Similarly, the current implementation splits the bounding box
   *    in half, a more optimal behaviour could be to split along the median
   *    of the data (as in quicksort etc).
   *
   * @tparam D the dimension to split on
   * @tparam DIM the number of dimensions
   * @tparam IndexIterator The data index iterator type
   * @tparam DataIterator The data iterator type
   * @tparam ListType The type of the collision pair list
   *
   * @param first The iterator pointing to the start of the range to partition
   * @param last The iterator pointing to the end of the range to partition
   * @param data The start of the data range
   * @param collisions The collision list
   * @param box The bouding box of the range we're to partition
   * @param depth The current recusion depth.
   */
  template <int D, int DIM, typename IndexIterator,
            typename DataIterator, typename ListType>
  void partition_data(IndexIterator first, IndexIterator last,
       DataIterator data, ListType &collisions,
       const bounding_box<DIM> &box, int depth)
  {
    // The next dimensions: X -> Y -> Z -> X ...
    const int D_NEXT = (D+1) % DIM;

    // Keep recusing until we either reach the maximum recusion depth or
    // the threshold of number of objects for brute force search is reached.
    if (depth + 1 < QC_MAXDEPTH && last - first > BF_THRESHOLD) {

      IndexIterator mid;

      // Copy the current level bounding box and starting with the lower
      // bound of the sub-division, split the box so that the max in the
      // current dimensio_along_axis is 1/2 of the current level box.
      bounding_box<DIM> sub_box(box);
      sub_box.max[D] = box.min[D] + (box.max[D] - box.min[D]) / 2;

      // Partition the space by the lower bound of the current dimension.
      // Then call the partition function on the next dimension, the order
      // of this is X -> Y -> Z -> X ...
      mid = std::partition(first, last, by_lower<D>(data, sub_box.max[D]));
      partition_data<D_NEXT>(first, mid, data, collisions, sub_box, depth+1);

      // Rejig the subdivided box to set the maximum of the current dimension
      // back to that of the current level box and the minimum to the maximum
      // of the lower bound box.
      sub_box.min[D] = sub_box.max[D];
      sub_box.max[D] = box.max[D];

      // Partition the space by the upper bound of the current dimension.
      // Then call the partition function on the next dimension, the order
      // of this is X -> Y -> Z -> X ...
      mid = std::partition(first, last, by_upper<D>(data, sub_box.min[D]));
      partition_data<D_NEXT>(mid, last, data, collisions, sub_box, depth+1);

    } else {

      // If the stopping condition has been met, then proceed with a
      // brute-force search for all the collison in the current level.
      detect_brute_force_w_check(first, last, data, collisions, box);
    }
  }

  /**
   * The "quick collide" aka "trevor" collision detection algorithm.
   *
   * The algorithm works by recursively splitting the space along each
   * dimension and partitioning the data that belongs in each side of the
   * split. Objects that span a split are handled by shuffling before
   * each side of the split if processed. Objects that collide are added
   * to a list of pairs of collisions.
   *
   * The list assumes there is a method "push_back" where new
   * elements can be added. Additionally it assumes that it defines a
   * dependant type "value_type" (which could be a std::pair<int, int>)
   * whose constructor takes two indices.
   *
   * @tparam Iterator The type of random access iterator
   * @tparam ListType a container that allows items to be appended.
   *
   * @param first The first iterator in the range
   * @param last The last iterator in the range
   * @param collisions The list of collisions
   */
  template <int DIM, typename Iterator, typename ListType>
  void detect_collisions(Iterator first, Iterator last, ListType &collisions) {

    // Ensure the amount of data is greater than zero.
    int n = last - first;
    assert(n > 0);

    // Create and fill a vector with the indices of the input data range
    std::vector<int> index(n);
    for (std::size_t i = 0; i < n; ++i) {
      index[i] = i;
    }

    // Get the bounding box of the whole data range
    bounding_box<DIM> box = get_bounding_box<DIM>(first, last);

    // Start the recursive partitioning of the data to find the collisions.
    partition_data<0>(index.begin(), index.end(), first, collisions, box, 0);
  }

  /** Wrapper function specialised for 2D collision detection */
  template <typename Iterator, typename ListType>
  void detect_collisions2d(Iterator first, Iterator last,
      ListType &collisions) {
    detect_collisions<2>(first, last, collisions);
  }

  /** Wrapper function specialised for 3D collision detection */
  template <typename Iterator, typename ListType>
  void detect_collisions3d(Iterator first, Iterator last,
      ListType &collisions) {
    detect_collisions<3>(first, last, collisions);
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPATIAL_INDEXING_DETECT_COLLISIONS_H
