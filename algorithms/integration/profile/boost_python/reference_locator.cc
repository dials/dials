/*
 * reference_locator.cc
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
#include <dials/algorithms/integration/profile/reference_locator.h>
#include <dials/algorithms/integration/profile/grid_sampler.h>
#include <dials/algorithms/integration/profile/xds_circle_sampler.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename ImageSampler>
  struct ReferenceLocatorPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const ReferenceLocator<ImageSampler> &r) {
      return boost::python::make_tuple(r.profile(), r.sampler());
    }
  };

  template <typename ImageSampler>
  void reference_locator_wrapper(const char *name) {
  
    typedef ImageSampler sampler_type;
    typedef ReferenceLocator<ImageSampler> locator_type;
  
    flex_double (locator_type::*profile_all)() const =
      &locator_type::profile;
    flex_double (locator_type::*profile_at_index)(
      std::size_t) const = &locator_type::profile;
    flex_double (locator_type::*profile_at_coord)(
      double3) const = &locator_type::profile;

    double3 (locator_type::*coord_at_index)(
      std::size_t) const = &locator_type::coord;
    double3 (locator_type::*coord_at_coord)(
      double3) const = &locator_type::coord;
  
    class_<locator_type>(name, no_init)
      .def(init<const flex_double &, const sampler_type&>((
        arg("profiles"), arg("sampler"))))
      .def("size", &locator_type::size)
      .def("sampler", &locator_type::sampler)
      .def("index", &locator_type::index)
      .def("profile", profile_all)
      .def("profile", profile_at_index)
      .def("profile", profile_at_coord)
      .def("coord", coord_at_index)
      .def("coord", coord_at_coord)
      .def("__len__", &locator_type::size)
      .def_pickle(ReferenceLocatorPickleSuite<sampler_type>());
  }

  void export_reference_locator()
  {
    reference_locator_wrapper<GridSampler>("GridReferenceLocator");
    reference_locator_wrapper<XdsCircleSampler>("XdsCircleReferenceLocator");
  }

}}} // namespace = dials::algorithms::boost_python
