/*
 * summation.cc
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
#include <dials/algorithms/integration/summation.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  /**
   * Integrate a reflection
   * @param r The reflection container
   */
  void summation3d(Reflection &r) {

    af::const_ref< int, af::c_grid<3> > shoebox_mask =
      r.get_shoebox_mask().const_ref();
    af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor());
    for (std::size_t i = 0; i < mask.size(); ++i) {
      mask[i] = (shoebox_mask[i] & shoebox::Valid) ? true : false;
    }

    // Integrate the reflection
    Summation result(r.get_shoebox().const_ref(),
                     r.get_shoebox_background().const_ref(),
                     mask.const_ref());

    // Set the intensity and variance
    r.set_intensity(result.intensity());
    r.set_intensity_variance(result.variance());
  }

  /**
   * Integrate a list of reflections
   * @param reflections The reflection list
   */
  void summation3d(af::ref<Reflection> reflections) {
    #pragma omp parallel for
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        if (reflections[i].is_valid()) {
          summation3d(reflections[i]);
        }
      } catch (dials::error) {
        reflections[i].set_valid(false);
      }
    }
  }

  void export_summation()
  {
    class_ <Summation> ("Summation", no_init)
      .def(init <const af::const_ref<double>&,
                 const af::const_ref<double>&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref<double>&,
                 const af::const_ref<double>&,
                 const af::const_ref<bool>&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< double, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< double, af::c_grid<2> >&,
                 const af::const_ref< bool, af::c_grid<2> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def(init <const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< double, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"))))
      .def(init <const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< double, af::c_grid<3> >&,
                 const af::const_ref< bool, af::c_grid<3> >&>((
          arg("signal"),
          arg("background"),
          arg("mask"))))
      .def("intensity", &Summation::intensity)
      .def("variance", &Summation::variance)
      .def("standard_deviation", &Summation::standard_deviation)
      .def("signal_intensity", &Summation::signal_intensity)
      .def("signal_variance", &Summation::signal_variance)
      .def("signal_standard_deviation", &Summation::signal_standard_deviation)
      .def("background_intensity", &Summation::background_intensity)
      .def("background_variance", &Summation::background_variance)
      .def("background_standard_deviation", 
        &Summation::background_standard_deviation);
        
    def("summation3d", (void(*)(Reflection&))&summation3d);
    def("summation3d", (void(*)(af::ref<Reflection>))&summation3d);
  }

}}} // namespace = dials::algorithms::boost_python
