/*
 * quadtree.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spatial_indexing/quadtree.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

//  template <typename ObjectType>
//  struct QuadTreeObject {
//    ObjectType object;
//    std::size_t index;
//    QuadTreeObject(const ObjectType &object_, std::size_t index_)
//      : object(object_),
//        index(index_) {}
//  };
//
//  template <>
//  bool contains<Box, QuadTreeObject>(const Box &a, const QuadTreeObject &b) {
//
//  }

//  template <>
//  bool collides<Box, QuadTreeObject>(const Box &a, const QuadTreeObject &b) {
//
//  }

//  template <typename ObjectType>
//  af::shared<ObjectType> quadtree_query_range(
//      const QuadTree<ObjectType> &tree,
//      const int4 &box) {
//    af::shared<ObjectType> result;
//    Box temp_box(box[0], box[2], box[1], box[3]);
//    bool success = tree.query_range(box, result);
//    DIALS_ASSERT(success);
//    return result;
//  }

//  template <typename ObjectType>
//  void quadtree_wrapper(const char *name) {

//    typedef ObjectType object_type;
//    typedef list_type::value_type value_type;
//    typedef QuadTree<value_type> quadtree_type;
//    typedef quadtree_type::box_type
//
//    class_<quadtree_type>(name, no_init)
//      def(init<quadtree_type::box_type,
//               quadtree_type::size_type>((
//        arg("box"),
//        arg("max_bucket_size")=10)))
//      : tree_(box, max_bucket_size)
//      .def("max_depth", &quadtree_type::max_depth)
//      .def("insert", &quadtree_type::insert)
//      .def("query_range", &quadtree_query_range)
//      .def("query_collision", &quadtree_query_collision);
//  }

  void export_quadtree()
  {
//    quadtree_wrapper< vec2<int> >("QuadTreeVec2int");
//    quadtree_wrapper< vec2<double> >("QuadTreeVec2Double");
  }

}}} // namespace = dials::algorithms::boost_python
