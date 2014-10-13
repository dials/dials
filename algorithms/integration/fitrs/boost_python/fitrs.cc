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
#include <dials/algorithms/integration/fitrs/fitrs.h>
#include <dials/algorithms/integration/fitrs/fit.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_profile_fitting_reciprocal_space()
  {
    vec3<double> (ProfileFittingReciprocalSpace::*call_single)(
      const TransformedShoebox&, vec3<double>) const =
        &ProfileFittingReciprocalSpace::operator();
    af::shared< vec3<double> > (ProfileFittingReciprocalSpace::*call_array)(
      const af::const_ref<TransformedShoebox>&,
      const af::const_ref< vec3<double> >&) const =
        &ProfileFittingReciprocalSpace::operator();

    class_<ProfileFittingReciprocalSpace>(
        "ProfileFittingReciprocalSpaceAlgorithm", no_init)
      .def(init<shared_ptr<ProfileFittingReciprocalSpace::locator_type>,
                std::size_t>((
        arg("locate"),
        arg("max_iter") = 10)))
      .def("__call__", call_single)
      .def("__call__", call_array);



    class_<Spec>("Spec", no_init)
      .def(init< const Beam&,
                 const Detector&,
                 const Goniometer&,
                 const Scan&,
                 double,
                 double >())
      ;


    class_<ReciprocalSpaceProfileFitting>("ReciprocalSpaceProfileFitting", no_init)
      .def(init< std::size_t >())
      .def("add", &ReciprocalSpaceProfileFitting::add)
      .def("execute", &ReciprocalSpaceProfileFitting::execute)
      ;
  }

}}} // namespace = dials::algorithms::boost_python
