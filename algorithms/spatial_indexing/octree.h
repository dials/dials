/*
 * octree.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPATIAL_INDEXING_OCTREE_H
#define DIALS_ALGORITHMS_SPATIAL_INDEXING_OCTREE_H

#include <cmath>
#include "qotree.h"
#include <dials/error.h>
#include <dials/algorithms/spatial_indexing/log2.h>

namespace dials { namespace algorithms {

  /** A 3d box structure */
  struct Box3d {
    int x0, y0, z0, x1, y1, z1;

    Box3d() : x0(0), y0(0), z0(0), x1(0), y1(0), z1(0) {}

    Box3d(int x0_, int y0_, int z0_, int x1_, int y1_, int z1_)
        : x0(x0_), y0(y0_), z0(z0_), x1(x1_), y1(y1_), z1(z1_) {}
  };

  /** Calculate the maximum depth */
  template <>
  std::size_t maximum_depth<Box3d>(const Box3d &box) {
    int w = box.x1 - box.x0;
    int h = box.y1 - box.y0;
    int d = box.z1 - box.z0;
    int min_d = (w < h ? (w < d ? w : d) : (h < d ? h : d));
    DIALS_ASSERT(min_d > 0);
    return (std::size_t)floor(log2(min_d));
  }

  /** Subdivide the box */
  template <>
  Box3d subdivide<Box3d>(const Box3d &box, std::size_t i) {
    int x0 = box.x0, x1 = box.x1, xc = (x0 + x1) / 2;
    int y0 = box.y0, y1 = box.y1, yc = (y0 + y1) / 2;
    int z0 = box.z0, z1 = box.z1, zc = (z0 + z1) / 2;
    int x00[8] = {x0, xc, x0, xc, x0, xc, x0, xc};
    int x11[8] = {xc, x1, xc, x1, xc, x1, xc, x1};
    int y00[8] = {y0, y0, yc, yc, y0, y0, yc, yc};
    int y11[8] = {yc, yc, y1, y1, yc, yc, y1, y1};
    int z00[8] = {z0, z0, z0, z0, zc, zc, zc, zc};
    int z11[8] = {zc, zc, zc, zc, z1, z1, z1, z1};
    return Box3d(x00[i], y00[i], z00[i], x11[i], y11[i], z11[i]);
  }

  /** Comparison operations for Box3d/Box3d objects */
  template <>
  struct compare<Box3d, Box3d> {
    /** Check if a box contains another box */
    static bool contains(const Box3d &box, const Box3d &v) {
      return box.x0 <= v.x0 && box.y0 <= v.y0 && box.z0 <= v.z0 && box.x1 >= v.x1
             && box.y1 >= v.y1 && box.z1 >= v.z1;
    }

    /** Check collisions between boxes */
    static bool collides(const Box3d &box, const Box3d &v) {
      return !(box.x0 >= v.x1 || v.x0 >= box.x1 || box.y0 >= v.y1 || v.y0 >= box.y1
               || box.z0 >= v.z1 || v.z0 >= box.z1);
    }
  };

  /**
   * A class implementing an octree. The class takes an object as template
   * parameter that is compared and stored in the octree.
   */
  template <typename ObjectType>
  class Octree {
  public:
    // General typedefs
    typedef QOTree<3, Box3d, ObjectType> tree_type;
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

    Octree(const box_type &box, size_type max_bucket_size = 10)
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

#endif  // DIALS_ALGORITHMS_SPATIAL_INDEXING_OCTREE_H
