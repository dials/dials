/*
 * experiment_ext.cc
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

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_beam();
  void export_goniometer();
  void export_detector();
  void export_scan();
  void export_detector_helpers();
  void export_scan_helpers();

  BOOST_PYTHON_MODULE(dials_model_experiment_ext)
  {
    export_beam();
    export_goniometer();
    export_detector();
    export_scan();
    export_detector_helpers();
    export_scan_helpers();
  }

}}} // namespace dials::model::boost_python
