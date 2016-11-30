/*
 * obs_manager_helpers.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <cctbx/miller.h>
#include <dials/algorithms/scaling/obs_manager_helpers.h>

namespace dials { namespace scaling { namespace boost_python {

  using namespace boost::python;

  typedef cctbx::miller::index<> miller_index;

  void export_minimum_multiplicity_selection()
  {
    def("minimum_multiplicity_selection", &minimum_multiplicity_selection, (
      arg("group_index"),
      arg("multiplicity")));
  }

  void export_grouped_obs()
  {
    // Create and return the wrapper for the grouped observations object
    class_ <GroupedObservations> ("GroupedObservations", no_init)
      .def(init< af::const_ref< miller_index >,
                 af::shared<double>,
                 af::shared<double>,
                 af::shared<double>,
                 af::shared<double> >((
        arg("group_index"),
        arg("intensity"),
        arg("weight"),
        arg("phi"),
        arg("scale"))))
      .def("get_intensity", &GroupedObservations::get_intensity)
      .def("set_intensity", &GroupedObservations::set_intensity,
        arg("intensity"))
      .def("get_weight", &GroupedObservations::get_weight)
      .def("set_weight", &GroupedObservations::set_weight,
        arg("weight"))
      .def("get_scale", &GroupedObservations::get_scale)
      .def("set_scale", &GroupedObservations::set_scale,
        arg("scale"))
      .def("get_group_index", &GroupedObservations::get_group_index)
      .def("get_phi", &GroupedObservations::get_phi)
      .def("get_group_size", &GroupedObservations::get_group_size)
      .def("get_average_intensity", &GroupedObservations::get_average_intensity)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

