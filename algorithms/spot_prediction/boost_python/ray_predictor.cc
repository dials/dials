/*
 * ray_predictor.cc
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
#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_ray_predictor()
  {
    // Useful typedefs
    typedef RayPredictor::reflection_type reflection_type;
    typedef RayPredictor::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (RayPredictor::*predict_single)(
      miller_index, mat3 <double>) const = &RayPredictor::operator();
    reflection_list_type (RayPredictor::*predict_array)(
      const flex_miller_index &, mat3 <double>) const = &RayPredictor::operator();

    // Create and return the wrapper for the spot predictor object
    class_ <RayPredictor> ("RayPredictor", no_init)
      .def(init <vec3<double>,
                 vec3<double>,
                 vec2<double> > ((
        arg("s0"),
        arg("m2"),
        arg("UB"),
        arg("dphi"))))
      .def("__call__", predict_single, (
        arg("miller_index"),
        arg("UB")))
      .def("__call__", predict_array, (
        arg("miller_indices"),
        arg("UB")));
  }

}}} // namespace = dials::spot_prediction::boost_python
