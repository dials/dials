/*
 * serialize_ext.cc
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

  void export_block_list();
  void export_serialize_shoebox();
  void export_shoebox_file_exporter();
  void export_shoebox_file_importer();

  BOOST_PYTHON_MODULE(dials_model_serialize_ext)
  {
    export_block_list();
    export_serialize_shoebox();
    export_shoebox_file_exporter();
    export_shoebox_file_importer();
  }

}}} // namespace = dials::model::boost_python
