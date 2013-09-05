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
#include <omptbx/omp_or_stubs.h>
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/algorithms/reflection_basis/rebin_pixels.h>
#include <dials/algorithms/reflection_basis/map_frames.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/algorithms/reflection_basis/map_pixels.h>
#include <dials/algorithms/reflection_basis/transform.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::int2;
  using scitbx::af::flex_grid;
  using scitbx::af::flex_int;
  using dials::model::ReflectionList;

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
    flex_vec3_double (*overload1)(const Detector&,
      const Beam&, std::size_t, bool) = &beam_vector_map;
    flex_vec3_double (*overload2)(const Detector&,
      const Beam&, bool) = &beam_vector_map;
    flex_vec3_double (*overload3)(const Detector&,
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

  void export_map_pixels()
  {
    class_<CoordinateGenerator>("CoordinateGenerator", no_init)
      .def(init<const CoordinateSystem&, int, int, const flex_vec3_double>((
          arg("cs"),
          arg("x0"),
          arg("y0"),
          arg("s1_map"))))
      .def("__call__", &CoordinateGenerator::operator());

    class_<GridIndexGenerator>("GridIndexGenerator", no_init)
      .def(init<const CoordinateSystem&, int, int, vec2<double>,
                std::size_t, const flex_vec3_double>((
          arg("cs"),
          arg("x0"),
          arg("y0"),
          arg("step_size"),
          arg("grid_half_size"),
          arg("s1_map"))))
      .def("__call__", &GridIndexGenerator::operator());

    flex_double (MapPixelsForward::*call_forward)(const CoordinateSystem&, int6,
        const flex_double&, const flex_bool&, const flex_double&) const =
          &MapPixelsForward::operator();

    flex_double (MapPixelsReverse::*call_reverse)(const CoordinateSystem&, int6,
        const flex_double&, const flex_double&) const =
          &MapPixelsReverse::operator();

    class_<MapPixelsForward>("MapPixelsForward", no_init)
      .def(init<const flex_vec3_double &,
                std::size_t,
                vec2<double> >((
        arg("s1_map"),
        arg("grid_half_size"),
        arg("step_size"))))
      .def("__call__", call_forward, (
        arg("cs"),
        arg("bbox"),
        arg("image"),
        arg("mask")));

    class_<MapPixelsReverse>("MapPixelsReverse", no_init)
      .def(init<const flex_vec3_double &,
                std::size_t,
                vec2<double> >((
        arg("s1_map"),
        arg("grid_half_size"),
        arg("step_size"))))
      .def("__call__", call_reverse, (
        arg("cs"),
        arg("bbox"),
        arg("grid")));
  }

  inline
  boost::shared_ptr<Forward> init_from_sweep_and_crystal(
      PyObject *sweep, PyObject *crystal,
      double nsigma, std::size_t grid_size) {
    return boost::shared_ptr<Forward>(new Forward(
      call_method<Beam>(sweep, "get_beam"),
      call_method<Detector>(sweep, "get_detector"),
      call_method<Goniometer>(sweep, "get_goniometer"),
      call_method<Scan>(sweep, "get_scan"),
      call_method<double>(crystal, "get_mosaicity", false),
      nsigma, grid_size));
  }

  inline
  void forward_with_reflection_list(const Forward &transform,
      ReflectionList &rlist) {
    // FIXME: Get Python error GC object already tracked. Possibly related to
    // creation of flex arrays in multiple threads
    //#pragma omp parallel for
    for (std::size_t i = 0; i < rlist.size(); ++i) {
      if (rlist[i].is_valid()) {
        try {
          flex_int shoebox_mask = rlist[i].get_shoebox_mask();
          flex_bool mask(shoebox_mask.accessor());
          for (std::size_t j = 0; j < mask.size(); ++j) {
            mask[j] = (shoebox_mask[j] & shoebox::Foreground);
          }
          flex_double grid = transform(
            rlist[i].get_beam_vector(),
            rlist[i].get_rotation_angle(),
            rlist[i].get_bounding_box(),
            rlist[i].get_shoebox(),
            mask);
          rlist[i].set_transformed_shoebox(grid);
//          std::pair<flex_double, flex_double> grid = transform(
//            rlist[i].get_beam_vector(),
//            rlist[i].get_rotation_angle(),
//            rlist[i].get_bounding_box(),
//            rlist[i].get_shoebox(),
//            rlist[i].get_shoebox_background(),
//            mask);
//          rlist[i].set_transformed_shoebox(grid.first);
////          rlist[i].set_transformed_shoebox_background(grid.second);
        } catch(dials::error) {
          rlist[i].set_valid(false);
        }
      }
    }
  }

  void export_transform()
  {
    flex_double(Forward::*forward_with_cs)(const CoordinateSystem&, int6,
      const flex_double&, const flex_bool&)const = &Forward::operator();
    flex_double(Forward::*forward_with_s1)(vec3<double>, double, int6,
      const flex_double&, const flex_bool&)const = &Forward::operator();
    std::pair<flex_double, flex_double>(Forward::*forward_with_background_cs)(
      const CoordinateSystem&, int6, const flex_double&, const flex_double&, 
      const flex_bool&)const = &Forward::operator();
    std::pair<flex_double, flex_double>(Forward::*forward_with_background_s1)(
      vec3<double>, double, int6, const flex_double&, const flex_double&, 
      const flex_bool&)const = &Forward::operator();

    flex_double(Reverse::*reverse_with_cs)(const CoordinateSystem&, int6,
      const flex_double&)const = &Reverse::operator();
    flex_double(Reverse::*reverse_with_s1)(vec3<double>, double, int6,
      const flex_double&)const = &Reverse::operator();

    boost_adaptbx::std_pair_conversions::to_tuple<flex_double, flex_double>();

    class_<Forward>("Forward", no_init)
      .def(init<const Beam&, const Detector&, const Goniometer&, const Scan&,
                double, double, std::size_t>((
        arg("beam"),
        arg("detector"),
        arg("gonio"),
        arg("scan"),
        arg("mosaicity"),
        arg("n_sigma"),
        arg("grid_half_size"))))
      .def("__init__",
        make_constructor(&init_from_sweep_and_crystal))
      .def("__call__", forward_with_cs, (
        arg("cs"),
        arg("bbox"),
        arg("image"),
        arg("mask")))
      .def("__call__", forward_with_s1, (
        arg("s1"),
        arg("phi"),
        arg("bbox"),
        arg("image"),
        arg("mask")))
      .def("__call__", forward_with_background_cs, (
        arg("cs"),
        arg("bbox"),
        arg("image"),
        arg("background"),
        arg("mask")))        
      .def("__call__", forward_with_background_s1, (
        arg("s1"),
        arg("phi"),
        arg("bbox"),
        arg("image"),
        arg("background"),
        arg("mask")))
      .def("__call__", &forward_with_reflection_list, (
        arg("rlist")));

    class_<Reverse>("Reverse", no_init)
      .def(init<const Beam&, const Detector&, const Goniometer&, const Scan&,
                double, double, std::size_t>((
        arg("beam"),
        arg("detector"),
        arg("gonio"),
        arg("scan"),
        arg("mosaicity"),
        arg("n_sigma"),
        arg("grid_half_size"))))
      .def("__call__",reverse_with_cs, (
        arg("cs"),
        arg("bbox"),
        arg("grid")))
      .def("__call__",reverse_with_s1, (
        arg("s1"),
        arg("phi"),
        arg("bbox"),
        arg("grid")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_reflection_basis_transform_ext)
  {
    export_rebin_pixels();
    export_map_frames();
    export_beam_vector_map();
    export_map_pixels();
    export_transform();
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
