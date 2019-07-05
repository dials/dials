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
#include <iterator>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <cmath>
#include <dials/error.h>
#include <dials/algorithms/spatial_indexing/log2.h>

namespace dials { namespace algorithms {

  using boost::mpl::for_each;
  using boost::mpl::range_c;

  //
  // User template specializations
  //

  /** Specify the data type of the bound coordinate */
  template <typename T>
  struct bound_coord_type {
    typedef double type;
  };

  /** Get the minimum bound of the box of type T in dimension DIM */
  template <int DIM, typename T>
  typename bound_coord_type<T>::type get_minimum_bound(const T &);

  /** Get the maximum bound of the box of type T in dimension DIM */
  template <int DIM, typename T>
  typename bound_coord_type<T>::type get_maximum_bound(const T &);

  //
  // BoundingBox struct and functions
  //

  /**
   * A simple bounding box structure that has min/max for each
   * dimension 0 <= D < DIM
   *
   * @tparam DIM the number of dimensions
   */
  template <int DIM, typename T>
  struct BoundingBox {
    static const int size = DIM;
    T min[DIM];
    T max[DIM];
    BoundingBox() {}
    BoundingBox(const BoundingBox &box) {
      for (std::size_t i = 0; i < DIM; ++i) {
        min[i] = box.min[i];
        max[i] = box.max[i];
      }
    }
  };

  /** The compile-time loop to set init the box elements */
  template <typename BoxType, typename ObjectType>
  struct init_bounding_box {
    const ObjectType &object_;
    BoxType &box_;

    init_bounding_box(const ObjectType &object, BoxType &box)
        : object_(object), box_(box) {}

    template <typename Index>
    void operator()(Index) {
      box_.min[Index::value] = get_minimum_bound<Index::value>(object_);
      box_.max[Index::value] = get_maximum_bound<Index::value>(object_);
    }
  };

  /** The compile-time loop to set min/max elements */
  template <typename BoxType, typename ObjectType>
  struct set_minmax_bounding_box {
    const ObjectType &object_;
    BoxType &box_;

    set_minmax_bounding_box(const ObjectType &object, BoxType &box)
        : object_(object), box_(box) {}

    template <typename Index>
    void operator()(Index) {
      if (get_minimum_bound<Index::value>(object_) < box_.min[Index::value]) {
        box_.min[Index::value] = get_minimum_bound<Index::value>(object_);
      }
      if (get_maximum_bound<Index::value>(object_) > box_.max[Index::value]) {
        box_.max[Index::value] = get_maximum_bound<Index::value>(object_);
      }
    }
  };

  /**
   * Get the bounding box from the input data. Basically scan the data for
   * the min/max in each dimension to get the bound that surrounds all the
   * elements. Use some template metaprogramming to cope with different
   * numbers of dimensions.
   */
  template <typename BoxType, typename Iterator>
  BoxType get_bounding_box(Iterator first, Iterator last) {
    typedef std::iterator_traits<Iterator> traits;
    typedef typename traits::value_type ObjectType;
    BoxType box;
    for_each<range_c<int, 0, BoxType::size> >(
      init_bounding_box<BoxType, ObjectType>(*first, box));
    for (Iterator it = first + 1; it < last; ++it) {
      for_each<range_c<int, 0, BoxType::size> >(
        set_minmax_bounding_box<BoxType, ObjectType>(*it, box));
    }
    return box;
  }

  //
  // BoxSize struct and functions
  //

  /** Hold a D-dimensional size */
  template <int DIM, typename T>
  struct BoxSize {
    static const int size = DIM;
    T d[DIM];
  };

  /** The compile-time loop to set init the size elements */
  template <typename SizeType, typename ObjectType>
  struct init_box_size {
    const ObjectType &object_;
    SizeType &size_;

    init_box_size(const ObjectType &object, SizeType &size)
        : object_(object), size_(size) {}

    template <typename Index>
    void operator()(Index) {
      size_.d[Index::value] = get_maximum_bound<Index::value>(object_)
                              - get_minimum_bound<Index::value>(object_);
    }
  };

  /** The compile-time loop to set min size elements */
  template <typename SizeType, typename ObjectType>
  struct set_minimum_box {
    const ObjectType &object_;
    SizeType &size_;

    set_minimum_box(const ObjectType &object, SizeType &size)
        : object_(object), size_(size) {}

    template <typename Index>
    void operator()(Index) {
      if (get_maximum_bound<Index::value>(object_)
            - get_minimum_bound<Index::value>(object_)
          < size_.d[Index::value]) {
        size_.d[Index::value] = get_maximum_bound<Index::value>(object_)
                                - get_minimum_bound<Index::value>(object_);
      }
    }
  };

  /** Loop through objects and find minimum size in each dimension */
  template <typename SizeType, typename Iterator>
  SizeType get_minimum_box_size(Iterator first, Iterator last) {
    typedef std::iterator_traits<Iterator> traits;
    typedef typename traits::value_type ObjectType;
    SizeType size;
    for_each<range_c<int, 0, SizeType::size> >(
      init_box_size<SizeType, ObjectType>(*first, size));
    for (Iterator it = first + 1; it < last; ++it) {
      for_each<range_c<int, 0, SizeType::size> >(
        set_minimum_box<SizeType, ObjectType>(*it, size));
    }
    return size;
  }

  //
  // DetectCollisions struct and functions
  //

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
  void detect_collisions_brute_force(Iterator first,
                                     Iterator last,
                                     ListType &collisions,
                                     Collides collides) {
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

  /** Template struct to check for collision */
  template <bool touching>
  struct no_collision_in_dimension;

  /**
   * Specialization of no collision struct. In this case boxes that touch
   * are considered to be colliding.
   */
  template <>
  struct no_collision_in_dimension<true> {
    template <typename CoordType>
    bool operator()(const CoordType &mina,
                    const CoordType &maxa,
                    const CoordType &minb,
                    const CoordType &maxb) {
      return mina > maxb || minb > maxa;
    }
  };

  /**
   * Specialization of no collision struct. In this case boxes that just touch
   * are not considered to be colliding.
   */
  template <>
  struct no_collision_in_dimension<false> {
    template <typename CoordType>
    bool operator()(const CoordType &mina,
                    const CoordType &maxa,
                    const CoordType &minb,
                    const CoordType &maxb) {
      return mina >= maxb || minb >= maxa;
    }
  };

  /** Compile time loop to iterate over dimensions to check for collisions */
  template <int I, bool touching>
  struct no_collision_loop {
    template <typename ObjectType>
    static bool exec(const ObjectType a, const ObjectType &b) {
      return no_collision_in_dimension<touching>()(get_minimum_bound<I - 1>(a),
                                                   get_maximum_bound<I - 1>(a),
                                                   get_minimum_bound<I - 1>(b),
                                                   get_maximum_bound<I - 1>(b))
             || no_collision_loop<I - 1, touching>::exec(a, b);
    }
  };

  /** Stop compile time loop to iterate over dimensions */
  template <bool touching>
  struct no_collision_loop<0, touching> {
    template <typename ObjectType>
    static bool exec(const ObjectType a, const ObjectType &b) {
      return false;
    }
  };

  /**
   * Check if the objects overlap (i.e. in 2D do the following)
   * return !(a.x0 > b.x1 || b.x0 > a.x1 || a.y0 > b.y1 || b.y0 > a.y1);
   *
   * To cope with different numbers of dimensions this is done using a bit
   * of template metaprogramming. The no_collision_loop struct sets up a
   * compile time loop over the number of dimensions.
   */
  template <int DIM, typename ObjectType, bool touching = true>
  struct collides {
    bool operator()(const ObjectType &a, const ObjectType &b) {
      return !no_collision_loop<DIM, touching>::exec(a, b);
    }
  };

  /**
   * Class to wrap up most of the collision detection code. Although the
   * user interface is through a function, the templatey nature of the
   * implementation means that it is a bit easier to understand what's going
   * on if we wrap it up in a class. The algorithm is first instantiated
   * and then called through the () operator.
   */
  template <int DIM, typename Iterator, typename ListType, bool touching = false>
  class DetectCollisions {
  public:
    // Types derived from data iterator
    typedef Iterator DataIterator;
    typedef std::iterator_traits<DataIterator> traits;
    typedef typename traits::value_type ObjectType;
    typedef typename bound_coord_type<ObjectType>::type CoordType;

    // Types derived from index list
    typedef std::vector<int> IndexList;
    typedef typename IndexList::value_type IndexType;
    typedef typename IndexList::iterator IndexIterator;

    // Types derived from collision list type
    typedef ListType CollisionList;
    typedef typename CollisionList::value_type Collision;

    // Bounding box and box size types
    typedef BoundingBox<DIM, CoordType> BoxType;
    typedef BoxSize<DIM, CoordType> DimType;

    // The collision checking object
    typedef collides<DIM, ObjectType, touching> collision_type;

    // The threshold to resort to a brute-force search
    static const int BF_THRESHOLD = 10;

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
     * @param first The first iterator in the range
     * @param last The last iterator in the range
     * @param collisions The list of collisions
     */
    void operator()(DataIterator first, DataIterator last, CollisionList &collisions) {
      // Ensure the amount of data is greater than zero.
      int n = last - first;
      DIALS_ASSERT(n > 0);

      // Create and fill a vector with the indices of the input data range
      IndexList index(n);
      for (std::size_t i = 0; i < n; ++i) {
        index[i] = i;
      }

      // Get the bounding box of the whole data range
      BoxType box = get_bounding_box<BoxType>(first, last);
      DimType min_size = get_minimum_box_size<DimType>(first, last);
      for (std::size_t i = 0; i < DIM; ++i) {
        DIALS_ASSERT(min_size.d[i] > 0);
      }

      // Calculate the maximum depth we can go to. Make sure that max depth
      // is atleast twice the size of the smallest object or zero.
      CoordType min_length = box.max[0] - box.min[0];
      std::size_t j = 0;
      for (std::size_t i = 0; i < DIM; ++i) {
        if (box.max[i] - box.min[i] < min_length) {
          min_length = box.max[i] - box.min[i];
          j = i;
        }
      }
      max_depth_ = log2(min_length / min_size.d[j]) - 1;
      if (max_depth_ < 1) max_depth_ = 1;
      max_depth_ *= DIM;

      // Start the recursive partitioning of the data to find the collisions.
      partition_data<0>(index.begin(), index.end(), first, collisions, box, 0);
    }

  private:
    int max_depth_;

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
     * half along the current axis. Parition the index array so that all
     * elements with the lower bound of their bounding box less than the split
     * are to the left and all those greater to the right. Then partition along
     * the next axis recursively.
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
     *    outliers. Similarly, the current implementation splits the bounding
     *    box in half, a more optimal behaviour could be to split along the
     *    median of the data (as in quicksort etc).
     *
     * @tparam D the dimension to split on
     *
     * @param first The iterator pointing to the start of the range to partition
     * @param last The iterator pointing to the end of the range to partition
     * @param data The start of the data range
     * @param collisions The collision list
     * @param box The bouding box of the range we're to partition
     * @param depth The current recusion depth.
     */
    template <int D>
    void partition_data(IndexIterator first,
                        IndexIterator last,
                        DataIterator data,
                        ListType &collisions,
                        const BoxType &box,
                        int depth) const {
      // The next dimensions: X -> Y -> Z -> X ...
      const int D_NEXT = (D + 1) % DIM;

      // Keep recusing until we either reach the maximum recusion depth or
      // the threshold of number of objects for brute force search is reached.
      if (depth < max_depth_ && last - first > BF_THRESHOLD) {
        IndexIterator mid;

        // Copy the current level bounding box and starting with the lower
        // bound of the sub-division, split the box so that the max in the
        // current dimensio_along_axis is 1/2 of the current level box.
        BoxType sub_box(box);
        sub_box.max[D] = box.min[D] + (box.max[D] - box.min[D]) / 2;

        // Partition the space by the lower bound of the current dimension.
        // Then call the partition function on the next dimension, the order
        // of this is X -> Y -> Z -> X ...
        mid = std::partition(first, last, by_lower<D>(data, sub_box.max[D]));
        partition_data<D_NEXT>(first, mid, data, collisions, sub_box, depth + 1);

        // Rejig the subdivided box to set the maximum of the current dimension
        // back to that of the current level box and the minimum to the maximum
        // of the lower bound box.
        sub_box.min[D] = sub_box.max[D];
        sub_box.max[D] = box.max[D];

        // Partition the space by the upper bound of the current dimension.
        // Then call the partition function on the next dimension, the order
        // of this is X -> Y -> Z -> X ...
        mid = std::partition(first, last, by_upper<D>(data, sub_box.min[D]));
        partition_data<D_NEXT>(mid, last, data, collisions, sub_box, depth + 1);

      } else {
        // If the stopping condition has been met, then proceed with a
        // brute-force search for all the collison in the current level.
        detect_brute_force_w_check(first, last, data, collisions, box);
      }
    }

    /**
     * Parition the data based on the lower bound of the data in the given
     * axis. All elements with their lower bound lower than the split to the
     * left, all others to the right.
     */
    template <int D>
    struct by_lower {
      const Iterator &data_;
      CoordType div_;

      by_lower(const DataIterator &data, CoordType div) : data_(data), div_(div) {}

      bool operator()(const IndexType &a) const {
        return get_minimum_bound<D>(*(data_ + a)) < div_;
      }
    };

    /**
     * Parition the data based on the upper bound of the data in the given
     * axis. All elements with their upper bound lower than the split to the
     * left, all others to the right.
     */
    template <int D>
    struct by_upper {
      const DataIterator &data_;
      CoordType div_;

      by_upper(const Iterator &data, CoordType div) : data_(data), div_(div) {}

      bool operator()(const int &a) const {
        return get_maximum_bound<D>(*(data_ + a)) < div_;
      }
    };

    /**
     * A helper function to wrap the call to detect_collisions_brute_force with
     * the check_collision struct.
     */
    void detect_brute_force_w_check(IndexIterator first,
                                    IndexIterator last,
                                    DataIterator data,
                                    ListType &collisions,
                                    const BoxType &box) const {
      // Do a brute force collision test using the collision checker to ensure
      // no pairs are added that have already been visited and added.
      detect_collisions_brute_force(
        first, last, collisions, check_collision(data, box));
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
    struct check_collision {
      const DataIterator d_;
      const BoxType &box_;

      /** Initialise the functor with the iterator and bounding box */
      check_collision(const DataIterator &d, const BoxType &box) : d_(d), box_(box) {}

      /**
       * Check for a collision. We return true if the objects collide and
       * if they have not already been set.
       */
      bool operator()(const IndexType &a, const IndexType &b) const {
        // If the objects collide and they have not already been
        // added then return true, otherwise return false.
        return collision_type()(*(d_ + a), *(d_ + b))
               && !already_visited(*(d_ + a), *(d_ + b), box_);
      }

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
      bool already_visited(const ObjectType &a,
                           const ObjectType &b,
                           const BoxType &box) const {
        bool result = false;
        for_each<range_c<int, 0, DIM> >(already_visited_loop(a, b, box, result));
        return result;
      }

      /** The compile-time loop to check if boxes have already been visited */
      struct already_visited_loop {
        const ObjectType &a_;
        const ObjectType &b_;
        const BoxType &box_;
        bool &result_;

        already_visited_loop(const ObjectType &a,
                             const ObjectType &b,
                             const BoxType &box,
                             bool &result)
            : a_(a), b_(b), box_(box), result_(result) {}

        template <typename I>
        void operator()(I) {
          result_ = result_
                    || (get_minimum_bound<I::value>(a_) < box_.min[I::value]
                        && get_minimum_bound<I::value>(b_) < box_.min[I::value]);
        }
      };
    };
  };

  /** Wrapper function specialised for 2D collision detection */
  template <typename Iterator, typename ListType>
  void detect_collisions2d(Iterator first, Iterator last, ListType &collisions) {
    // Create the collision detection object and call
    DetectCollisions<2, Iterator, ListType, false>()(first, last, collisions);
  }

  /** Wrapper function specialised for 3D collision detection */
  template <typename Iterator, typename ListType>
  void detect_collisions3d(Iterator first, Iterator last, ListType &collisions) {
    // Create the collision detection object and call
    DetectCollisions<3, Iterator, ListType, false>()(first, last, collisions);
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPATIAL_INDEXING_DETECT_COLLISIONS_H
