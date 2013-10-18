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
#include <dials/algorithms/reflection_basis/map_frames.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/algorithms/reflection_basis/index_generator.h>
#include <dials/algorithms/reflection_basis/transform.h>
#include <dials/algorithms/reflection_basis/ideal_profile.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;
  using dials::model::Reflection;
  using scitbx::af::int2;
  using scitbx::af::flex_double;

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

  inline
  boost::shared_ptr<TransformSpec> init_forward_from_sweep_and_crystal(
      PyObject *sweep, PyObject *crystal,
      double nsigma, std::size_t grid_size) {
    return boost::shared_ptr<TransformSpec>(new TransformSpec(
      call_method<Beam>(sweep, "get_beam"),
      call_method<Detector>(sweep, "get_detector"),
      call_method<Goniometer>(sweep, "get_goniometer"),
      call_method<Scan>(sweep, "get_scan"),
      call_method<double>(crystal, "get_mosaicity", false),
      nsigma, grid_size));
  }

  void export_forward() {

    class_<TransformSpec>("TransformSpec", no_init)
      .def(init<const Beam &, const Detector &, 
                const Goniometer &, const Scan &, 
                double, double, std::size_t>())
      .def("__init__", make_constructor(&init_forward_from_sweep_and_crystal))
      .def("m2", &TransformSpec::m2)
      .def("s0", &TransformSpec::s0)
      .def("image_size", &TransformSpec::image_size)
      .def("grid_size", &TransformSpec::grid_size)
      .def("step_size", &TransformSpec::step_size)
      .def("grid_centre", &TransformSpec::grid_centre)
      .def("s1_map", &TransformSpec::s1_map);

    class_<Forward>("Forward", no_init)
      .def(init<const TransformSpec&,
                const vec3<double>&, double, int6,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec&,
                const vec3<double>&, double, int6,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec&,
                const CoordinateSystem&, int6,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec&,
                const CoordinateSystem&, int6,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< double, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec&,
                const Reflection&,
                bool>())
      .def("profile", &Forward::profile)
      .def("background", &Forward::background)
      .def("zfraction", &Forward::zfraction);
      
    def("forward_batch", &forward_batch);
  }

  void export_ideal_profile() {
    def("ideal_profile", &ideal_profile);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_reflection_basis_transform_ext)
  {
    export_map_frames();
    export_beam_vector_map();
    export_index_generator();
    export_forward();
    export_ideal_profile();
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
