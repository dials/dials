/*
 * model_ext.cc
 *
 *   Copyright (C) 2013 Diamond Light Source, James Parkhurst
 *
 *   This code is distributed under the BSD license, a copy of which is
 *   included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials { namespace model { namespace experiment { namespace boost_python {

  void export_beam();
  void export_goniometer();
  //void export_detector();

  BOOST_PYTHON_MODULE(dials_model_experiment_ext)
  {
    export_beam();
    export_goniometer();
    //export_detector();
  }

}}}} // namespace dials::model::experiment::boost_python
