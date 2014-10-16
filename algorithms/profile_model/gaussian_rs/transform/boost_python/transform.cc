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
#include <dials/algorithms/profile_model/gaussian_rs/transform.h>

namespace dials { namespace algorithms { namespace profile_model { namespace gaussian_rs { {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  boost::shared_ptr< TransformSpec<FloatType> >
    init_forward_from_experiment(
      object experiment, double sigma_b, double sigma_m,
      double nsigma, std::size_t grid_size) {
    return boost::shared_ptr< TransformSpec<FloatType> >(
      new TransformSpec<FloatType>(
        extract<Beam>(experiment.attr("beam")),
        extract<Detector>(experiment.attr("detector")),
        extract<Goniometer>(experiment.attr("goniometer")),
        extract<Scan>(experiment.attr("scan")),
        sigma_b, sigma_m, nsigma, grid_size));
  }

  template <typename FloatType>
  void forward_wrapper(const char *name) {

    typedef Forward<FloatType> ForwardType;
    typedef TransformSpec<FloatType> TransformSpecType;

    class_<TransformSpecType>("TransformSpec", no_init)
      .def(init<const Beam &, const Detector &,
                const Goniometer &, const Scan &,
                double, double, double, std::size_t>())
      .def("__init__", make_constructor(
        &init_forward_from_experiment<FloatType>))
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
    forward_wrapper<double>("Forward");
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
