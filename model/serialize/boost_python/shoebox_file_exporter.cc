/*
 * shoebox_file_exporter.cc
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
#include <dials/model/serialize/shoebox_file_exporter.h>

namespace dials { namespace model { namespace serialize {
  namespace boost_python {

  using namespace boost::python;

  void export_shoebox_file_exporter()
  {
    class_<ShoeboxFileExporter, boost::noncopyable>(
        "ShoeboxFileExporter", no_init)
      .def(init<const std::string &,
                const af::const_ref<std::size_t>&,
                const af::const_ref<int6>,
                const af::const_ref<double>,
                std::size_t,
                std::size_t,
                const std::string&>((
        arg("filename"),
        arg("panel"),
        arg("bbox"),
        arg("num_frame"),
        arg("num_panel"),
        arg("blob") = "")))
      .def("flush", &ShoeboxFileExporter::flush)
      .def("finished", &ShoeboxFileExporter::finished)
      .def("next", &ShoeboxFileExporter::next);
  }

}}}} // namespace dials::model::serialize::boost_python

