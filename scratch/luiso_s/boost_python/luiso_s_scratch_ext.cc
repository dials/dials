/*
 * luiso_s_scratch_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: Luiso & James
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace scratch { namespace boost_python {
  using namespace boost::python;

  void luiso_s_scratch_ext();

  BOOST_PYTHON_MODULE(luiso_s_scratch_ext)
  {
    luiso_s_scratch_ext();
  }

}}}
