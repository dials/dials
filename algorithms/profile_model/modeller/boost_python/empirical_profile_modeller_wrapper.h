

#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_BOOST_PYTHON_EMPIRICAL_PROFILE_MODELLER_WRAPPER_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_BOOST_PYTHON_EMPIRICAL_PROFILE_MODELLER_WRAPPER_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename T>
  class_<T, bases<ProfileModellerIface> > empirical_profile_modeller_wrapper(
    const char *name) {
    class_<T, bases<ProfileModellerIface> > result(name, no_init);
    /* bases<ProfileModellerIface> >("EmpiricalProfileModeller", no_init) */
    result.def("add", &T::add)
      .def("valid", &T::valid)
      .def("n_reflections", &T::n_reflections)
      .def("model", &T::model)
      .def("fit", &T::fit)
      .def("validate", &T::validate)
      .def("accumulate", &T::accumulate)
      //.def("finalize", &T::finalize)
      .def("finalized", &T::finalized)
      .def("data", &T::data)
      .def("mask", &T::mask)
      .def("size", &T::size)
      .def("copy", &T::copy)
      .def("__len__", &T::size);

    return result;
  }

}}}  // namespace dials::algorithms::boost_python

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_BOOST_PYTHON_EMPIRICAL_PROFILE_MODELLER_WRAPPER_H
