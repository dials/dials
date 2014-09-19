/*
 * fit_ext.cc
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
#include <dials/algorithms/integration/fit_image/fit.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_integration_fit_image_ext)
  {
    typedef ImageSpaceProfileFitting::reference_learner_type reference_learner_type;

    class_<Spec>("Spec", no_init)
      .def(init< const Beam&,
                 const Detector&,
                 const Goniometer&,
                 const Scan&,
                 double,
                 double >())
      ;

    class_<reference_learner_type>("ReferenceLearner", no_init)
      .def("get", &reference_learner_type::get)
      .def("data", &reference_learner_type::data)
      .def("count", &reference_learner_type::count)
      .def("__len__", &reference_learner_type::size)
      ;

    class_<ImageSpaceProfileFitting>("ImageSpaceProfileFitting", no_init)
      .def(init< std::size_t >())
      .def("add", &ImageSpaceProfileFitting::add)
      .def("execute", &ImageSpaceProfileFitting::execute)
      ;
  }

}}} // namespace = dials::algorithms::boost_python
