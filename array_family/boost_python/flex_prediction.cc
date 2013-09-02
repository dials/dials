/*
 * flex_prediction.cc
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
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/model/data/prediction.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using dials::model::Prediction;

  void export_flex_prediction()
  {
    scitbx::af::boost_python::flex_wrapper <
        Prediction, return_internal_reference<> >::plain("prediction");
  }

}}} // namespace dials::af::boost_python
