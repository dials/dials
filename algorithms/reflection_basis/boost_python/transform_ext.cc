/*
 * transform_ext.cc
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

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  void export_map_frames();
  void export_beam_vector_map();
  void export_index_generator();
  void export_ideal_profile();
  void export_transform();

  BOOST_PYTHON_MODULE(dials_algorithms_reflection_basis_transform_ext)
  {
    export_map_frames();
    export_beam_vector_map();
    export_index_generator();
    export_ideal_profile();
    export_transform();
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
