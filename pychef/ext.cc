/*
 * ext.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <dials/pychef/Chef.h>

namespace dials { namespace pychef { namespace boost_python {

  using namespace boost::python;

  void export_observations() {
    class_<Observations::map_type>("ObservationGroupMap")
      .def(map_indexing_suite<Observations::map_type>());

    class_<Observations>("Observations", no_init)
      .def(init<scitbx::af::const_ref<cctbx::miller::index<> > const &,
                sgtbx::space_group,
                bool>((arg("miller_index"), arg("space_group"), arg("anomalous_flag"))))
      .def("observation_groups", &Observations::observation_groups);
  }

  void export_observation_group() {
    class_<ObservationGroup>("ObservationGroup", no_init)
      .def(init<cctbx::miller::index<>, bool>((arg("miller_index"), arg("flag"))))
      .def("add_iplus", &ObservationGroup::add_iplus)
      .def("add_iminus", &ObservationGroup::add_iminus)
      .def("miller_index", &ObservationGroup::miller_index)
      .def("iplus", &ObservationGroup::iplus)
      .def("iminus", &ObservationGroup::iminus)
      .def("is_centric", &ObservationGroup::is_centric);
  }

  void export_chef_statistics() {
    typedef ChefStatistics chef_statistics_t;
    class_<chef_statistics_t>("ChefStatistics", no_init)
      .def(init<scitbx::af::const_ref<cctbx::miller::index<> > const &,
                af::const_ref<double> const &,
                af::const_ref<double> const &,
                af::const_ref<double> const &,
                af::const_ref<std::size_t> const &,
                af::const_ref<std::size_t> const &,
                cctbx::miller::binner const &,
                sgtbx::space_group,
                bool,
                int>((arg("miller_index"),
                      arg("intensities"),
                      arg("sigmas"),
                      arg("d_star_sq"),
                      arg("dose"),
                      arg("counts_complete"),
                      arg("binner"),
                      arg("space_group"),
                      arg("anomalous_flag"),
                      arg("n_steps"))))
      .def("iplus_completeness", &chef_statistics_t::iplus_completeness)
      .def("iminus_completeness", &chef_statistics_t::iminus_completeness)
      .def("ieither_completeness", &chef_statistics_t::ieither_completeness)
      .def("iboth_completeness", &chef_statistics_t::iboth_completeness)
      .def("iplus_completeness_bins", &chef_statistics_t::iplus_completeness_bins)
      .def("iminus_completeness_bins", &chef_statistics_t::iminus_completeness_bins)
      .def("ieither_completeness_bins", &chef_statistics_t::ieither_completeness_bins)
      .def("iboth_completeness_bins", &chef_statistics_t::iboth_completeness_bins)
      .def("rcp_bins", &chef_statistics_t::rcp_bins)
      .def("scp_bins", &chef_statistics_t::scp_bins)
      .def("rcp", &chef_statistics_t::rcp)
      .def("scp", &chef_statistics_t::scp)
      .def("rd", &chef_statistics_t::rd);
  }

  BOOST_PYTHON_MODULE(dials_pychef_ext) {
    export_observations();
    export_observation_group();
    export_chef_statistics();
  }

}}}  // namespace dials::pychef::boost_python
