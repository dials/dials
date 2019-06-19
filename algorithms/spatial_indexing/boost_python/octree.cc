/*
 * octree.cc
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
#include <scitbx/vec3.h>
#include <dials/algorithms/spatial_indexing/octree.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  namespace boost_python {

    template <typename ObjectType>
    struct OctreeObject {
      ObjectType object;
      std::size_t index;
      OctreeObject(const ObjectType &object_, std::size_t index_)
          : object(object_), index(index_) {}
    };

  }  // namespace boost_python

  template <typename Object>
  struct compare<Box3d, boost_python::OctreeObject<Object> > {
    static bool contains(const Box3d &box,
                         const boost_python::OctreeObject<Object> &v) {
      return box.x0 <= v.object[0] && box.y0 <= v.object[1] && box.z0 <= v.object[2]
             && box.x1 >= v.object[0] && box.y1 >= v.object[1] && box.z1 >= v.object[2];
    }

    static bool collides(const Box3d &box,
                         const boost_python::OctreeObject<Object> &v) {
      return contains(box, v);
    }
  };

  namespace boost_python {

    using namespace boost::python;
    using scitbx::vec3;

    template <typename ObjectType>
    Octree<OctreeObject<ObjectType> > *make_octree(af::int6 box, int max_bucket_size) {
      return new Octree<OctreeObject<ObjectType> >(
        Box3d(box[0], box[2], box[4], box[1], box[3], box[5]), max_bucket_size);
    }

    template <typename ObjectType>
    af::shared<std::size_t> octree_query_range(
      const Octree<OctreeObject<ObjectType> > &tree,
      const af::int6 &box) {
      af::shared<std::size_t> result;
      af::shared<OctreeObject<ObjectType> > temp;
      bool success =
        tree.query_range(Box3d(box[0], box[2], box[4], box[1], box[3], box[5]), temp);
      DIALS_ASSERT(success);
      for (std::size_t i = 0; i < temp.size(); ++i) {
        result.push_back(temp[i].index);
      }
      return result;
    }

    template <typename ObjectType>
    af::shared<std::size_t> octree_query_collision(
      const Octree<OctreeObject<ObjectType> > &tree,
      const ObjectType &obj) {
      af::shared<std::size_t> result;
      af::shared<OctreeObject<ObjectType> > temp;
      OctreeObject<ObjectType> obj2(obj, 0);
      bool success = tree.query_collision(obj2, temp);
      DIALS_ASSERT(success);
      for (std::size_t i = 0; i < temp.size(); ++i) {
        result.push_back(temp[i].index);
      }
      return result;
    }

    template <typename ObjectType>
    Octree<OctreeObject<ObjectType> >
      *make_octree_from_list(af::const_ref<ObjectType> a, int max_bucket_size) {
      DIALS_ASSERT(a.size() > 0);
      typename ObjectType::value_type xmin, xmax, ymin, ymax, zmin, zmax;
      xmin = xmax = a[0][0];
      ymin = ymax = a[0][1];
      zmin = zmax = a[0][2];
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i][0] < xmin) xmin = a[i][0];
        if (a[i][0] > xmax) xmax = a[i][0];
        if (a[i][1] < ymin) ymin = a[i][1];
        if (a[i][1] > ymax) ymax = a[i][1];
        if (a[i][2] < zmin) zmin = a[i][2];
        if (a[i][2] > zmax) zmax = a[i][2];
      }
      af::int6 box(
        (int)xmin, (int)xmax + 1, (int)ymin, (int)ymax + 1, (int)zmin, (int)zmax + 1);
      Octree<OctreeObject<ObjectType> > *qt =
        make_octree<ObjectType>(box, max_bucket_size);
      for (std::size_t i = 0; i < a.size(); ++i) {
        qt->insert(OctreeObject<ObjectType>(a[i], i));
      }
      return qt;
    }

    template <typename ObjectType>
    void octree_wrapper(const char *name) {
      typedef OctreeObject<ObjectType> object_type;
      typedef Octree<object_type> octree_type;

      class_<octree_type>(name, no_init)
        .def("__init__",
             make_constructor(&make_octree_from_list<ObjectType>,
                              default_call_policies(),
                              (arg("items"), arg("max_bucket_size") = 10)))
        .def("max_depth", &octree_type::max_depth)
        .def("query_range",
             &octree_query_range<ObjectType>);  //.def("query_collision",
                                                //&octree_query_collision<ObjectType>);

      def("make_spatial_index",
          &make_octree_from_list<ObjectType>,
          return_value_policy<manage_new_object>(),
          (arg("items"), arg("max_bucket_size") = 10));
    }

    void export_octree() {
      octree_wrapper<vec3<int> >("OctreeVec3int");
      octree_wrapper<vec3<double> >("OctreeVec3Double");
    }

  }  // namespace boost_python
}}   // namespace dials::algorithms
