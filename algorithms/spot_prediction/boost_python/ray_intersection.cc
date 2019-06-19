/*
 * ray_intersection.cc
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
#include <dxtbx/model/detector.h>
#include <dials/algorithms/spot_prediction/ray_intersection.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_ray_intersection() {
    af::shared<bool> (*ray_intersection_table)(const Detector&, af::reflection_table) =
      &ray_intersection;
    af::shared<bool> (*ray_intersection_table_w_panel)(
      const Detector&, af::reflection_table, std::size_t panel) = &ray_intersection;
    af::shared<bool> (*ray_intersection_table_w_panels)(
      const Detector&, af::reflection_table, af::ref<std::size_t> panel) =
      &ray_intersection;

    // init<
    // const Detector&,
    // const af::const_ref< vec3<double> >&
    //> from_s1_any_panel((arg("detector"), arg("s1")));

    // init<
    // const Detector&,
    // const af::const_ref< vec3<double> >&,
    // std::size_t
    //> from_s1_single_panel((arg("detector"), arg("s1"), arg("panel")));

    // init<
    // const Detector&,
    // const af::const_ref< vec3<double> >&,
    // const af::const_ref<std::size_t>&
    //> from_s1_panel_array((arg("detector"), arg("s1"), arg("panel")));

    // class_<ray_intersection2>("ray_intersection2", no_init)
    //.def(from_s1_any_panel)
    //.def(from_s1_single_panel)
    //.def(from_s1_panel_array);

    // Export all the ray intersection functions
    def("ray_intersection",
        ray_intersection_table,
        (arg("detector"), arg("reflection_table")));
    def("ray_intersection",
        ray_intersection_table_w_panel,
        (arg("detector"), arg("reflection_table"), arg("panel")));
    def("ray_intersection",
        ray_intersection_table_w_panels,
        (arg("detector"), arg("reflection_table"), arg("panel")));
  }

}}}  // namespace dials::algorithms::boost_python
