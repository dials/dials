/*
 * indexing_ext.cc
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fft3d();

  BOOST_PYTHON_MODULE(dials_algorithms_indexing_ext)
  {
    export_fft3d();
  }

}}} // namespace = dials::algorithms::boost_python
