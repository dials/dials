/*
 * beam_vector_map.cc
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
#include <dials/algorithms/reflection_basis/beam_vector_map.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  void export_beam_vector_map()
  {
    af::versa< vec3<double>, af::c_grid<2> > (*overload1)(const Detector&,
      const Beam&, std::size_t, bool) = &beam_vector_map;
    af::versa< vec3<double>, af::c_grid<2> > (*overload2)(const Detector&,
      const Beam&, bool) = &beam_vector_map;
    af::versa< vec3<double>, af::c_grid<2> > (*overload3)(const Detector&,
      const Beam&) = &beam_vector_map;

    def("beam_vector_map", overload1, (
      arg("detector"),
      arg("beam"),
      arg("n_div"),
      arg("corners")));
    def("beam_vector_map", overload2, (
      arg("detector"),
      arg("beam"),
      arg("corners")));
    def("beam_vector_map", overload3, (
      arg("detector"),
      arg("beam")));
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
