/*
 * quadtree.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPATIAL_INDEXING_QUADTREE_H
#define DIALS_ALGORITHMS_SPATIAL_INDEXING_QUADTREE_H

#include <cmath>
#include "qotree.h"
#include <dials/error.h>
#include <dials/algorithms/spatial_indexing/log2.h>

namespace dials { namespace algorithms {

  /** A 2d box structure */
  struct Box {
    int x0, y0, x1, y1;

    Box() : x0(0), y0(0), x1(0), y1(0) {}

    Box(int x0_, int y0_, int x1_, int y1_) : x0(x0_), y0(y0_), x1(x1_), y1(y1_) {}
  };

  /** Calculate the maximum depth */
  template <>
  std::size_t maximum_depth<Box>(const Box &box) {
    int w = box.x1 - box.x0;
    int h = box.y1 - box.y0;
    int min_d = w < h ? w : h;
    DIALS_ASSERT(min_d > 0);
    return (std::size_t)floor(log2(min_d));
  }

  /** Subdivide the box */
  template <>
  Box subdivide<Box>(const Box &box, std::size_t i) {
    int x0 = box.x0, x1 = box.x1, xc = (x0 + x1) / 2;
    int y0 = box.y0, y1 = box.y1, yc = (y0 + y1) / 2;
    int x00[4] = {x0, xc, x0, xc};
    int x11[4] = {xc, x1, xc, x1};
    int y00[4] = {y0, y0, yc, yc};
    int y11[4] = {yc, yc, y1, y1};
    return Box(x00[i], y00[i], x11[i], y11[i]);
  }

  /** Comparison operations for Box/Box objects */
  template <>
  struct compare<Box, Box> {
    /** Check if a box contains another box */
    static bool contains(const Box &box, const Box &v) {
      return box.x0 <= v.x0 && box.y0 <= v.y0 && box.x1 >= v.x1 && box.y1 >= v.y1;
    }

    /** Check collisions between boxes */
    static bool collides(const Box &box, const Box &v) {
      return !(box.x0 >= v.x1 || v.x0 >= box.x1 || box.y0 >= v.y1 || v.y0 >= box.y1);
    }
  };

  /**
   * A class implementing a quadtree. The class takes an object as template
   * parameter that is compared and stored in the quadtree.
   */
  template <typename ObjectType>
  class Quadtree {
  public:
    // General typedefs
    typedef QOTree<2, Box, ObjectType> tree_type;
    typedef typename tree_type::box_type box_type;
    typedef typename tree_type::object_type object_type;
    typedef typename tree_type::size_type size_type;

    // Object list typedefs
    typedef typename tree_type::object_list_type object_list_type;
    typedef typename tree_type::object_iterator object_iterator;
    typedef typename tree_type::const_object_iterator const_object_iterator;

    // Node list typedefs
    typedef typename tree_type::node_type node_type;
    typedef typename tree_type::node_list_type node_list_type;
    typedef typename tree_type::node_pointer node_pointer;
    typedef typename tree_type::const_node_pointer const_node_pointer;
    typedef typename tree_type::node_iterator node_iterator;
    typedef typename tree_type::const_node_iterator const_node_iterator;

    Quadtree(const box_type &box, size_type max_bucket_size = 10)
        : tree_(box, max_bucket_size) {}

    /** Get the const iterator to the beginning of the node list */
    const_node_iterator node_begin() const {
      return tree_.node_begin();
    }

    /** Get the const iterator to the end of the node list */
    const_node_iterator node_end() const {
      return tree_.node_end();
    }

    /** Get the maximum bucket size */
    size_type max_bucket_size() const {
      return tree_.max_bucket_size();
    }

    /** Get the maximum depth of the tree */
    size_type max_depth() const {
      return tree_.max_depth();
    }

    /** Insert an element into the tree */
    bool insert(const object_type &v) {
      return tree_.insert(v);
    }

    /** Find all the objects in the given range */
    template <typename ObjectList>
    bool query_range(const box_type &range, ObjectList &elements) const {
      return tree_.query_range(range, elements);
    }

    /** Find all the objects colliding with the given object */
    template <typename ObjectList>
    bool query_collision(const object_type &v, ObjectList &elements) const {
      return tree_.query_collision(v, elements);
    }

  private:
    tree_type tree_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPATIAL_INDEXING_QUADTREE_H
