/*
 * xds_circle_sampler.cc
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
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  class XdsCircleSamplerIterator {
  public:
    XdsCircleSamplerIterator(const XdsCircleSampler &obj, int index)
      : obj_(obj), index_(index) {}
      
    typedef double2 value_type;
    typedef int difference_type;
    typedef double2 reference;
    typedef const double2* pointer;
    typedef std::input_iterator_tag iterator_category;
    value_type operator*() { return obj_[index_]; }
    XdsCircleSamplerIterator operator++() { ++index_; return *this; }
    XdsCircleSamplerIterator operator--() { --index_; return *this; }
    XdsCircleSamplerIterator operator++(int) { 
      XdsCircleSamplerIterator result(obj_, index_);
      ++(*this); 
      return result; 
    }
    XdsCircleSamplerIterator operator--(int) {
      XdsCircleSamplerIterator result(obj_, index_);
      --(*this); 
      return result; 
    }
    bool operator==(const XdsCircleSamplerIterator& rhs) { 
      return index_ == rhs.index_; 
    }
    bool operator!=(const XdsCircleSamplerIterator& rhs) { 
      return index_ == rhs.index_; 
    }
  private:
    const XdsCircleSampler &obj_;
    int index_;
  };
  
  XdsCircleSamplerIterator xds_sampler_begin(const XdsCircleSampler &obj) {
    return XdsCircleSamplerIterator(obj, 0);
  }

  XdsCircleSamplerIterator xds_sampler_end(const XdsCircleSampler &obj) {
    return XdsCircleSamplerIterator(obj, obj.size());
  }

  void export_xds_circle_sampler()
  {
    class_<XdsCircleSampler>("XdsCircleSampler", no_init)
      .def(init<int2>((
        arg("image_size"))))
      .def("image_size", &XdsCircleSampler::image_size)
      .def("image_centre", &XdsCircleSampler::image_centre)
      .def("r0", &XdsCircleSampler::r0)
      .def("r1", &XdsCircleSampler::r1)
      .def("r2", &XdsCircleSampler::r2)
      .def("nearest", &XdsCircleSampler::nearest)
      .def("size", &XdsCircleSampler::size)
      .def("__getitem__", &XdsCircleSampler::operator[])
      .def("__len__", &XdsCircleSampler::size)
      .def("__iter__", range(xds_sampler_begin, xds_sampler_end));
  }

}}} // namespace = dials::algorithms::boost_python
