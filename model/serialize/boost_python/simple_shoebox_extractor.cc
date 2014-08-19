/*
 * simple_shoebox_extractor.cc
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
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/model/serialize/simple_shoebox_extractor.h>

namespace dials { namespace model { namespace serialize {
  namespace boost_python {

  using namespace boost::python;

  void export_simple_shoebox_extractor()
  {
    class_<SimpleShoeboxExtractor>("SimpleShoeboxExtractor", no_init)
      .def(init<
          const af::const_ref< Shoebox<> >&,
          int,
          int,
          std::size_t>())
      .def("next", &SimpleShoeboxExtractor::next)
      ;
  }

}}}} // namespace dials::model::serialize::boost_python

