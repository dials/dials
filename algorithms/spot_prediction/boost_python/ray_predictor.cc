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
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static
  af::reflection_table call_with_miller_index_array(
      const ScanStaticRayPredictor &self,
      const af::const_ref< cctbx::miller::index<> > &h,
      const mat3<double> &UB) {
    af::reflection_table result;
    af::shared< cctbx::miller::index<> > hkl = result["miller_index"];
    af::shared< vec3<double> > s1 = result["s1"];
    af::shared< bool > entering = result["entering"];
    af::shared< double > phi = result["phi"];
    for (std::size_t i = 0; i < h.size(); ++i) {
      af::small<Ray,2> rays = self(h[i], UB);
      for (std::size_t j = 0; j < rays.size(); ++j) {
        hkl.push_back(h[i]);
        s1.push_back(rays[i].s1);
        entering.push_back(rays[i].entering);
        phi.push_back(rays[i].angle);
      }
    }
    return result;
  }

  void export_ray_predictor()
  {
    // Useful typedefs
    typedef RayPredictor::reflection_type reflection_type;
    typedef RayPredictor::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (RayPredictor::*predict_single)(
      miller_index, mat3 <double>) const = &RayPredictor::operator();
    reflection_list_type (RayPredictor::*predict_array)(
      const af::const_ref<miller_index> &, mat3 <double>) const =
        &RayPredictor::operator();

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

    // Create and return the wrapper for the spot predictor object
    class_ <ScanStaticRayPredictor> ("ScanStaticRayPredictor", no_init)
      .def(init <vec3<double>,
                 vec3<double>,
                 vec2<double> > ((
        arg("s0"),
        arg("m2"),
        arg("dphi"))))
      .def("__call__", &ScanStaticRayPredictor::operator(), (
        arg("miller_index"),
        arg("UB")))
      .def("__call__", &call_with_miller_index_array);
  }

}}} // namespace = dials::spot_prediction::boost_python
