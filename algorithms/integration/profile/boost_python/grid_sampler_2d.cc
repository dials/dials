/*
 * grid_sampler_2d.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/profile/grid_sampler_2d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  class GridSampler2DIterator {
  public:
    GridSampler2DIterator(const GridSampler2D &obj, int index)
      : obj_(obj), index_(index) {}

    typedef double2 value_type;
    typedef int difference_type;
    typedef double2 reference;
    typedef const double2* pointer;
    typedef std::input_iterator_tag iterator_category;
    value_type operator*() { return obj_[index_]; }
    GridSampler2DIterator operator++() { ++index_; return *this; }
    GridSampler2DIterator operator--() { --index_; return *this; }
    GridSampler2DIterator operator++(int) {
      GridSampler2DIterator result(obj_, index_);
      ++(*this);
      return result;
    }
    GridSampler2DIterator operator--(int) {
      GridSampler2DIterator result(obj_, index_);
      ++(*this);
      return result;
    }
    bool operator==(const GridSampler2DIterator& rhs) {
      return index_ == rhs.index_;
    }
    bool operator!=(const GridSampler2DIterator& rhs) {
      return index_ == rhs.index_;
    }
  private:
    const GridSampler2D &obj_;
    int index_;
  };

  GridSampler2DIterator grid_sampler_begin(const GridSampler2D &obj) {
    return GridSampler2DIterator(obj, 0);
  }

  GridSampler2DIterator grid_sampler_end(const GridSampler2D &obj) {
    return GridSampler2DIterator(obj, obj.size());
  }


  struct GridSampler2DPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const GridSampler2D &obj) {
      return boost::python::make_tuple(obj.image_size(), obj.grid_size());
    }
  };

  void export_grid_sampler_2d()
  {
    class_<GridSampler2D>("GridSampler2D", no_init)
      .def(init<int2, int2>((
        arg("image_size"),
        arg("grid_size"))))
      .def("image_size", &GridSampler2D::image_size)
      .def("grid_size", &GridSampler2D::grid_size)
      .def("step_size", &GridSampler2D::step_size)
      .def("nearest", &GridSampler2D::nearest)
      .def("nearest_n", &GridSampler2D::nearest_n)
      .def("size", &GridSampler2D::size)
      .def("weight", &GridSampler2D::weight)
      .def("neighbours", &GridSampler2D::neighbours)
      .def("__getitem__", &GridSampler2D::operator[])
      .def("__len__", &GridSampler2D::size)
      .def("__iter__", range(&grid_sampler_begin, &grid_sampler_end))
      .def_pickle(GridSampler2DPickleSuite());
  }

}}} // namespace = dials::algorithms::boost_python

