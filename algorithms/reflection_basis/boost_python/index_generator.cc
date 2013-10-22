/*
 * index_generator.cc
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
#include <dials/algorithms/reflection_basis/index_generator.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  CoordinateGenerator* make_coordinate_generator(
      const CoordinateSystem &cs, 
      int x0, int y0, 
      af::versa< vec3<double>, scitbx::af::flex_grid<> > s1_map) {
    return new CoordinateGenerator(cs, x0, y0, 
      af::versa< vec3<double>, af::c_grid<2> >(s1_map.handle(), 
        af::c_grid<2>(s1_map.accessor()))); 
  }

  GridIndexGenerator* make_grid_index_generator(const CoordinateSystem &cs, 
      int x0, int y0, vec2<double> step_size, std::size_t grid_half_size, 
      af::versa< vec3<double>, scitbx::af::flex_grid<> > s1_map) {
    return new GridIndexGenerator(cs, x0, y0, step_size, grid_half_size, 
      af::versa< vec3<double>, af::c_grid<2> >(
        s1_map.handle(), af::c_grid<2>(s1_map.accessor()))); 
  }

  void export_index_generator()
  {
    class_<CoordinateGenerator>("CoordinateGenerator", no_init)
      .def("__init__", make_constructor(
        &make_coordinate_generator, 
        default_call_policies(), (
          arg("cs"),
          arg("x0"),
          arg("y0"),
          arg("s1_map"))))
      .def("__call__", &CoordinateGenerator::operator());

    class_<GridIndexGenerator>("GridIndexGenerator", no_init)
      .def("__init__", make_constructor(
        &make_grid_index_generator, 
        default_call_policies(), (
          arg("cs"),
          arg("x0"),
          arg("y0"),
          arg("step_size"),
          arg("grid_half_size"),
          arg("s1_map"))))    
      .def("__call__", &GridIndexGenerator::operator());
  }
  
}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
