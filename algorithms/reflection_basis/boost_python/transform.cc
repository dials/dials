/*
 * transform.cc
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
#include <dials/algorithms/reflection_basis/transform.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  boost::shared_ptr< TransformSpec<FloatType> >
    init_forward_from_sweep_and_crystal(
      PyObject *sweep, PyObject *crystal,
      double nsigma, std::size_t grid_size) {
    return boost::shared_ptr< TransformSpec<FloatType> >(
      new TransformSpec<FloatType>(
        call_method<Beam>(sweep, "get_beam"),
        call_method<Detector>(sweep, "get_detector"),
        call_method<Goniometer>(sweep, "get_goniometer"),
        call_method<Scan>(sweep, "get_scan"),
        call_method<double>(crystal, "get_mosaicity", false),
        nsigma, grid_size));
  }

  template <typename FloatType>
  void forward_wrapper(const char *name) {

    typedef Forward<FloatType> ForwardType;
    typedef TransformSpec<FloatType> TransformSpecType;

    class_<TransformSpecType>("TransformSpec", no_init)
      .def(init<const Beam &, const Detector &,
                const Goniometer &, const Scan &,
                double, double, std::size_t>())
      .def("__init__", make_constructor(
        &init_forward_from_sweep_and_crystal<FloatType>))
      .def("m2", &TransformSpecType::m2)
      .def("s0", &TransformSpecType::s0)
      .def("image_size", &TransformSpecType::image_size)
      .def("grid_size", &TransformSpecType::grid_size)
      .def("step_size", &TransformSpecType::step_size)
      .def("grid_centre", &TransformSpecType::grid_centre)
      .def("s1_map", &TransformSpecType::s1_map);

    class_<ForwardType>(name, no_init)
      .def(init<const TransformSpec<FloatType>&,
                const vec3<double>&, double, int6,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec<FloatType>&,
                const vec3<double>&, double, int6,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec<FloatType>&,
                const CoordinateSystem&, int6,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def(init<const TransformSpec<FloatType>&,
                const CoordinateSystem&, int6,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< FloatType, af::c_grid<3> >&,
                const af::const_ref< bool, af::c_grid<3> >& >())
      .def("profile", &ForwardType::profile)
      .def("background", &ForwardType::background)
      .def("zfraction", &ForwardType::zfraction);
  }

  void export_transform() {
    forward_wrapper<Shoebox<>::float_type>("Forward");
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
