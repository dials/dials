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
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/vec2.h>
#include <dials/algorithms/spatial_indexing/quadtree.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  namespace boost_python {

    template <typename ObjectType>
    struct QuadtreeObject {
      ObjectType object;
      std::size_t index;
      QuadtreeObject(const ObjectType &object_, std::size_t index_)
          : object(object_), index(index_) {}
    };

  }  // namespace boost_python

  template <typename Object>
  struct compare<Box, boost_python::QuadtreeObject<Object> > {
    static bool contains(const Box &box,
                         const boost_python::QuadtreeObject<Object> &v) {
      return box.x0 <= v.object[0] && box.y0 <= v.object[1] && box.x1 >= v.object[0]
             && box.y1 >= v.object[1];
    }

    static bool collides(const Box &box,
                         const boost_python::QuadtreeObject<Object> &v) {
      return contains(box, v);
    }
  };

  namespace boost_python {

    using namespace boost::python;
    using scitbx::vec2;

    template <typename ObjectType>
    Quadtree<QuadtreeObject<ObjectType> > *make_quadtree(af::int4 box,
                                                         int max_bucket_size) {
      return new Quadtree<QuadtreeObject<ObjectType> >(
        Box(box[0], box[2], box[1], box[3]), max_bucket_size);
    }

    template <typename ObjectType>
    af::shared<std::size_t> quadtree_query_range(
      const Quadtree<QuadtreeObject<ObjectType> > &tree,
      const af::int4 &box) {
      af::shared<std::size_t> result;
      af::shared<QuadtreeObject<ObjectType> > temp;
      bool success = tree.query_range(Box(box[0], box[2], box[1], box[3]), temp);
      DIALS_ASSERT(success);
      for (std::size_t i = 0; i < temp.size(); ++i) {
        result.push_back(temp[i].index);
      }
      return result;
    }

    template <typename ObjectType>
    af::shared<std::size_t> quadtree_query_collision(
      const Quadtree<QuadtreeObject<ObjectType> > &tree,
      const ObjectType &obj) {
      af::shared<std::size_t> result;
      af::shared<QuadtreeObject<ObjectType> > temp;
      QuadtreeObject<ObjectType> obj2(obj, 0);
      bool success = tree.query_collision(obj2, temp);
      DIALS_ASSERT(success);
      for (std::size_t i = 0; i < temp.size(); ++i) {
        result.push_back(temp[i].index);
      }
      return result;
    }

    template <typename ObjectType>
    Quadtree<QuadtreeObject<ObjectType> >
      *make_quadtree_from_list(af::const_ref<ObjectType> a, int max_bucket_size) {
      DIALS_ASSERT(a.size() > 0);
      typename ObjectType::value_type xmin, xmax, ymin, ymax;
      xmin = xmax = a[0][0];
      ymin = ymax = a[0][1];
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i][0] < xmin) xmin = a[i][0];
        if (a[i][0] > xmax) xmax = a[i][0];
        if (a[i][1] < ymin) ymin = a[i][1];
        if (a[i][1] > ymax) ymax = a[i][1];
      }
      af::int4 box((int)xmin, (int)xmax + 1, (int)ymin, (int)ymax + 1);
      Quadtree<QuadtreeObject<ObjectType> > *qt =
        make_quadtree<ObjectType>(box, max_bucket_size);
      for (std::size_t i = 0; i < a.size(); ++i) {
        qt->insert(QuadtreeObject<ObjectType>(a[i], i));
      }
      return qt;
    }

    template <typename ObjectType>
    void quadtree_wrapper(const char *name) {
      typedef QuadtreeObject<ObjectType> object_type;
      typedef Quadtree<object_type> quadtree_type;

      class_<quadtree_type>(name, no_init)
        .def("__init__",
             make_constructor(&make_quadtree_from_list<ObjectType>,
                              default_call_policies(),
                              (arg("items"), arg("max_bucket_size") = 10)))
        .def("max_depth", &quadtree_type::max_depth)
        .def(
          "query_range",
          &quadtree_query_range<ObjectType>);  //.def("query_collision",
                                               //&quadtree_query_collision<ObjectType>);

      def("make_spatial_index",
          &make_quadtree_from_list<ObjectType>,
          return_value_policy<manage_new_object>(),
          (arg("items"), arg("max_bucket_size") = 10));
    }

    void export_quadtree() {
      quadtree_wrapper<vec2<int> >("QuadtreeVec2int");
      quadtree_wrapper<vec2<double> >("QuadtreeVec2Double");
    }

  }  // namespace boost_python
}}   // namespace dials::algorithms
