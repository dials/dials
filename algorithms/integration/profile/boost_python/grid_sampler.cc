/*
 * grid_sampler.cc
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
#include <dials/algorithms/integration/profile/grid_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  class GridSamplerIterator {
  public:
    GridSamplerIterator(const GridSampler &obj, int index)
      : obj_(obj), index_(index) {}
      
    typedef double3 value_type;
    typedef int difference_type;
    typedef double3 reference;
    typedef const double3* pointer;
    typedef std::input_iterator_tag iterator_category;
    value_type operator*() { return obj_[index_]; }
    GridSamplerIterator operator++() { ++index_; return *this; }
    GridSamplerIterator operator--() { --index_; return *this; }
    GridSamplerIterator operator++(int) {
      GridSamplerIterator result(obj_, index_);
      ++(*this); 
      return result; 
    }
    GridSamplerIterator operator--(int) {
      GridSamplerIterator result(obj_, index_);
      ++(*this); 
      return result; 
    }
    bool operator==(const GridSamplerIterator& rhs) { 
      return index_ == rhs.index_; 
    }
    bool operator!=(const GridSamplerIterator& rhs) { 
      return index_ == rhs.index_; 
    }
  private:
    const GridSampler &obj_;
    int index_;
  };
  
  GridSamplerIterator grid_sampler_begin(const GridSampler &obj) {
    return GridSamplerIterator(obj, 0);
  }

  GridSamplerIterator grid_sampler_end(const GridSampler &obj) {
    return GridSamplerIterator(obj, obj.size());
  }


  struct GridSamplerPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const GridSampler &obj) {
      return boost::python::make_tuple(obj.volume_size(), obj.grid_size());
    }
  };

  void export_grid_sampler()
  {
    class_<GridSampler>("GridSampler", no_init)
      .def(init<int3, int3>((
        arg("volume_size"),
        arg("grid_size"))))
      .def("image_size", &GridSampler::volume_size)
      .def("grid_size", &GridSampler::grid_size)
      .def("step_size", &GridSampler::step_size)
      .def("nearest", &GridSampler::nearest)
      .def("size", &GridSampler::size)
      .def("__getitem__", &GridSampler::operator[])
      .def("__len__", &GridSampler::size)
      .def("__iter__", range(&grid_sampler_begin, &grid_sampler_end))
      .def_pickle(GridSamplerPickleSuite());
  }

}}} // namespace = dials::algorithms::boost_python
