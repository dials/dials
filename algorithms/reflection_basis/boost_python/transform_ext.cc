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
#include <scitbx/array_family/flex_types.h>
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
  using dials::model::Reflection;
  using scitbx::af::int2;
  using scitbx::af::flex_double;

  inline
  af::versa< double, af::c_grid<2> > rebin_pixels_wrapper(
      const af::const_ref< double, af::c_grid<2> > &input,
      const af::const_ref< vec2<double>, af::c_grid<2> > &inputxy, int2 size) {
    af::versa< double, af::c_grid<2> > output(af::c_grid<2>(size[0], size[1]), 
      af::init_functor_null<double>());
    af::ref< double, af::c_grid<2> > output_ref = output.ref();
    rebin_pixels(output_ref, input, inputxy);
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

  MapPixelsForward* make_map_pixels_forward(
      af::versa< vec3<double>, scitbx::af::flex_grid<> > s1_map,
      std::size_t grid_half_size, vec2<double> step_size) {
    return new MapPixelsForward(af::versa< vec3<double>, af::c_grid<2> >(
        s1_map.handle(), af::c_grid<2>(s1_map.accessor())),
        grid_half_size, step_size); 
  }

  MapPixelsReverse* make_map_pixels_reverse(
      af::versa< vec3<double>, scitbx::af::flex_grid<> > s1_map,
      std::size_t grid_half_size, vec2<double> step_size) {
    return new MapPixelsReverse(af::versa< vec3<double>, af::c_grid<2> >(
        s1_map.handle(), af::c_grid<2>(s1_map.accessor())),
        grid_half_size, step_size); 
  }

  void export_map_pixels()
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

    af::versa< double, af::c_grid<3> > (MapPixelsForward::*call_forward)(
        const CoordinateSystem&, int6, 
        const af::const_ref< double, af::c_grid<3> >&, 
        const af::const_ref< bool, af::c_grid<3> >&, 
        const af::const_ref< double, af::c_grid<2> >&) const = 
          &MapPixelsForward::operator();

    af::versa< double, af::c_grid<3> > (MapPixelsReverse::*call_reverse)(
        const CoordinateSystem&, int6,
        const af::const_ref< double, af::c_grid<3> >&, 
        const af::const_ref< double, af::c_grid<2> >&) const =
          &MapPixelsReverse::operator();

    class_<MapPixelsForward>("MapPixelsForward", no_init)
      .def("__init__", make_constructor(
        &make_map_pixels_forward, 
        default_call_policies(), (
          arg("s1_map"),  
          arg("grid_half_size"),
          arg("step_size"))))
      .def("__call__", call_forward, (
        arg("cs"),
        arg("bbox"),
        arg("image"),
        arg("mask")));

    class_<MapPixelsReverse>("MapPixelsReverse", no_init)
      .def("__init__", make_constructor(
        &make_map_pixels_reverse, 
        default_call_policies(), (
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
      af::ref<Reflection> rlist) {
    // FIXME: Get Python error GC object already tracked. Possibly related to
    // creation of shared arrays in multiple threads
    //#pragma omp parallel for
    for (std::size_t i = 0; i < rlist.size(); ++i) {
      if (rlist[i].is_valid()) {
        try {
          af::versa< int, af::c_grid<3> > shoebox_mask = rlist[i].get_shoebox_mask();
          af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor(), 
            af::init_functor_null<bool>());
          for (std::size_t j = 0; j < mask.size(); ++j) {
            mask[j] = (shoebox_mask[j] & shoebox::Foreground);
          }
          af::versa< double, af::c_grid<3> > grid = transform(
            rlist[i].get_beam_vector(),
            rlist[i].get_rotation_angle(),
            rlist[i].get_bounding_box(),
            rlist[i].get_shoebox().const_ref(),
            mask.const_ref());
          rlist[i].set_transformed_shoebox(grid);
//          std::pair<versa< double, c_grid<3> >, versa< double, c_grid<3> > > 
//            grid = transform(
//              rlist[i].get_beam_vector(),
//              rlist[i].get_rotation_angle(),
//              rlist[i].get_bounding_box(),
//              rlist[i].get_shoebox(),
//              rlist[i].get_shoebox_background(),
//              mask);
//          rlist[i].set_transformed_shoebox(grid.first);
////          rlist[i].set_transformed_shoebox_background(grid.second);
        } catch(dials::error) {
          rlist[i].set_valid(false);
        }
      }
    }
  }
  
  std::pair<flex_double, flex_double> forward_with_background_cs(
      const Forward &transform,
      const CoordinateSystem& cs, int6 bbox, 
      const af::const_ref< double, af::c_grid<3> > &shoebox, 
      const af::const_ref< double, af::c_grid<3> > &background, 
      const af::const_ref< bool, af::c_grid<3> > &mask) {
  
    std::pair< af::versa< double, af::c_grid<3> >, 
               af::versa< double, af::c_grid<3> > >
      result = transform(cs, bbox, shoebox, background, mask);
      return std::make_pair(flex_double(
        result.first.handle(), result.first.accessor().as_flex_grid()),
        flex_double(result.second.handle(), 
          result.second.accessor().as_flex_grid()));
  }
  
  std::pair<flex_double, flex_double> forward_with_background_s1(
      const Forward &transform,
      vec3<double> s1, double phi, int6 bbox, 
      const af::const_ref< double, af::c_grid<3> > &shoebox, 
      const af::const_ref< double, af::c_grid<3> > &background, 
      const af::const_ref< bool, af::c_grid<3> > &mask) {
  
    std::pair< af::versa< double, af::c_grid<3> >, 
               af::versa< double, af::c_grid<3> > >
      result = transform(s1, phi, bbox, shoebox, background, mask);
      return std::make_pair(flex_double(
        result.first.handle(), result.first.accessor().as_flex_grid()),
        flex_double(result.second.handle(), 
          result.second.accessor().as_flex_grid()));
  }

  void export_transform()
  {
    af::versa< double, af::c_grid<3> >(Forward::*forward_with_cs)(
      const CoordinateSystem&, int6,
      const af::const_ref< double, af::c_grid<3> >&, 
      const af::const_ref< bool, af::c_grid<3> >&)const = &Forward::operator();
    af::versa< double, af::c_grid<3> >(Forward::*forward_with_s1)(
      vec3<double>, double, int6,
      const af::const_ref< double, af::c_grid<3> >&, 
      const af::const_ref< bool, af::c_grid<3> >&)const = &Forward::operator();
    af::versa< double, af::c_grid<3> >(Reverse::*reverse_with_cs)(
      const CoordinateSystem&, int6,
      const af::const_ref< double, af::c_grid<3> >&)const = &Reverse::operator();
    af::versa< double, af::c_grid<3> >(Reverse::*reverse_with_s1)(
      vec3<double>, double, int6,
      const af::const_ref< double, af::c_grid<3> >&)const = &Reverse::operator();

    boost_adaptbx::std_pair_conversions::to_tuple<
      scitbx::af::flex_double, scitbx::af::flex_double>();

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
