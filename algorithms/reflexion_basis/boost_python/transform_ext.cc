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
#include <dials/algorithms/reflexion_basis/rebin_pixels.h>
#include <dials/algorithms/reflexion_basis/map_frames.h>
#include <dials/algorithms/reflexion_basis/beam_vector_map.h>

namespace dials { namespace algorithms { namespace reflexion_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::int2;
  using scitbx::af::flex_grid;
  
  inline
  flex_double rebin_pixels_wrapper(const flex_double &input, 
      const flex_vec2_double &inputxy, int2 size) {
    flex_double output(flex_grid<>(size[0], size[1]));
    rebin_pixels(output, input, inputxy);
    return output;
  }

  void export_rebin_pixels() 
  {
    def("rebin_pixels", &rebin_pixels_wrapper, (
      arg("input"), arg("inputxy"), arg("size")));
  }
  
  void export_map_frames()
  {
    class_<MapFramesForward>(
        "MapFramesForward", no_init)
      .def(init<double, double, double, double, int>((
          arg("starting_angle"),
          arg("oscillation"),
          arg("mosaicity"),
          arg("n_sigma"),
          arg("grid_size_e3"))))
      .def("__call__",
        &MapFramesForward::operator(), (
          arg("frames"),
          arg("phi"),
          arg("zeta")));

    class_<MapFramesReverse>(
        "MapFramesReverse", no_init)
      .def(init<double, double, double, double, int>((
          arg("starting_angle"),
          arg("oscillation"),
          arg("mosaicity"),
          arg("n_sigma"),
          arg("grid_size_e3"))))
      .def("__call__",
        &MapFramesReverse::operator(), (
          arg("frames"),
          arg("phi"),
          arg("zeta")));
  }
  
  void export_beam_vector_map()
  {
    flex_vec3_double (*overload1)(const Detector&, const Beam&, std::size_t, bool) = &beam_vector_map;
    flex_vec3_double (*overload2)(const Detector&, const Beam&, bool) = &beam_vector_map;
    flex_vec3_double (*overload3)(const Detector&, const Beam&) = &beam_vector_map;
  
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

  BOOST_PYTHON_MODULE(dials_algorithms_reflexion_basis_transform_ext)
  {
    export_rebin_pixels();
    export_map_frames();
  }

}}}}} // namespace = dials::algorithms::reflexion_basis::transform::boost_python
