/*
 * shoebox_file_importer.cc
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
#include <dials/model/serialize/shoebox_file_importer.h>

namespace dials { namespace model { namespace serialize {
  namespace boost_python {

  using namespace boost::python;

  ShoeboxFileImporter* make_shoebox_file_importer(
      const std::string &filename,
      boost::python::tuple gmt,
      boost::python::tuple dmt,
      boost::python::tuple mmt) {
 
    // The input
    ShoeboxFileImporter::gain_map_array_type gm;
    ShoeboxFileImporter::dark_map_array_type dm;
    ShoeboxFileImporter::mask_map_array_type mm;

    // Ensure all the same length
    std::size_t num = len(gmt);
    DIALS_ASSERT(len(dmt) == num);
    DIALS_ASSERT(len(mmt) == num);

    // Extract the stuff
    for (std::size_t i = 0; i < num; ++i) {
      gm.push_back(extract<ShoeboxFileImporter::gain_map_type>(gmt[i]));
      dm.push_back(extract<ShoeboxFileImporter::dark_map_type>(dmt[i]));
      mm.push_back(extract<ShoeboxFileImporter::mask_map_type>(mmt[i]));
    }

    // Return the new importer
    return new ShoeboxFileImporter(filename, 
      gm.const_ref(), dm.const_ref(), mm.const_ref());
  }

  void export_shoebox_file_importer()
  {
    class_<ShoeboxFileImporter, boost::noncopyable>(
        "ShoeboxFileImporter", no_init)
      .def("init", make_constructor(
        &make_shoebox_file_importer,
        default_call_policies(), (
          arg("filename"),
          arg("gain"),
          arg("dark"),
          arg("mask"))))
      ;
  }

}}}} // namespace dials::model::serialize::boost_python


