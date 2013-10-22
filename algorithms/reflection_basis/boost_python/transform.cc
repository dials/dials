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
  
  static
  Forward<Reflection::float_type> forward_with_reflection(
      const TransformSpec<Reflection::float_type> &spec,
      const Reflection &reflection, bool with_background) {

    // Get reflection properties
    vec3<double> s1 = reflection.get_beam_vector();
    double phi = reflection.get_rotation_angle();
    int6 bbox = reflection.get_bounding_box();
    
    // Get the shoebox data
    af::const_ref< Reflection::float_type, af::c_grid<3> > shoebox = 
      reflection.get_shoebox().const_ref();
    af::const_ref< Reflection::float_type, af::c_grid<3> > background = 
      reflection.get_shoebox_background().const_ref();
    af::const_ref< int, af::c_grid<3> > shoebox_mask =
      reflection.get_shoebox_mask().const_ref();

    // Create the mask
    af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor(), false);
    for (std::size_t i = 0; i < mask.size(); ++i) {
      mask[i] = shoebox_mask[i] & Valid && shoebox_mask[i] & Foreground;
    }

    // Do the transform
    if (with_background) {
      return Forward<Reflection::float_type>(spec, s1, phi, bbox, shoebox, 
        background, mask.const_ref());
    } else {
      return Forward<Reflection::float_type>(spec, s1, phi, bbox, shoebox, 
        mask.const_ref());
    }
  }

  static
  void forward_batch(const TransformSpec<Reflection::float_type> &spec,
                     af::ref<Reflection> rlist, bool with_background) {
    for (std::size_t i = 0; i < rlist.size(); ++i) {
      if (rlist[i].is_valid()) {
        Forward<Reflection::float_type> transform = forward_with_reflection(
          spec, rlist[i], with_background);
        rlist[i].set_transformed_shoebox(transform.profile());
        if (with_background) {
          rlist[i].set_transformed_shoebox_background(transform.background());
        }
      }
    }
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

    forward_wrapper<Reflection::float_type>("Forward");
      
    def("forward", &forward_with_reflection);
    def("forward_batch", &forward_batch);
  }
  
}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python
