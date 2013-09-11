/*
 * profile_fitting_reciprocal_space.cc
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
#include <dials/algorithms/integration/profile_fitting_reciprocal_space.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_profile_fitting_reciprocal_space()
  {
    void (ProfileFittingReciprocalSpace::*call_single)(
      Reflection&) const = &ProfileFittingReciprocalSpace::operator();
    void (ProfileFittingReciprocalSpace::*call_array)(
      af::ref<Reflection>) const = &ProfileFittingReciprocalSpace::operator();
  
    class_<ProfileFittingReciprocalSpace>(
        "ProfileFittingReciprocalSpaceAlgorithm", no_init)
      .def(init<shared_ptr<ProfileFittingReciprocalSpace::locator_type> >((
        arg("locate"))))
      .def("__call__", call_single)
      .def("__call__", call_array);
  }

}}} // namespace = dials::algorithms::boost_python
