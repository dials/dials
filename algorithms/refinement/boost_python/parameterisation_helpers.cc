
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include "../parameterisation/parameterisation_helpers.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_parameterisation_helpers() {
    def("multi_panel_compose",
        &multi_panel_compose,
        (arg("initial_state"),
         arg("params_vals"),
         arg("params_axes"),
         arg("detector"),
         arg("offsets"),
         arg("dir1s"),
         arg("dir2s"),
         arg("Tau1"),
         arg("dTau1_dtau1"),
         arg("Tau1"),
         arg("dTau1_dtau1"),
         arg("Tau1"),
         arg("dTau1_dtau1")));

    def("multi_panel_compose",
        &selected_multi_panel_compose,
        (arg("initial_state"),
         arg("params_vals"),
         arg("params_axes"),
         arg("detector"),
         arg("selection"),
         arg("offsets"),
         arg("dir1s"),
         arg("dir2s"),
         arg("Tau1"),
         arg("dTau1_dtau1"),
         arg("Tau1"),
         arg("dTau1_dtau1"),
         arg("Tau1"),
         arg("dTau1_dtau1")));

    class_<CrystalOrientationCompose>("CrystalOrientationCompose", no_init)
      .def(init<mat3<double>,
                double,
                vec3<double>,
                double,
                vec3<double>,
                double,
                vec3<double> >((arg("U0"),
                                arg("phi1"),
                                arg("phi1_axis"),
                                arg("phi2"),
                                arg("phi2_axis"),
                                arg("phi3"),
                                arg("phi3_axis"))))
      .def("U", &CrystalOrientationCompose::U)
      .def("dU_dphi1", &CrystalOrientationCompose::dU_dphi1)
      .def("dU_dphi2", &CrystalOrientationCompose::dU_dphi2)
      .def("dU_dphi3", &CrystalOrientationCompose::dU_dphi3);

    class_<PanelGroupCompose>("PanelGroupCompose", no_init)
      .def(init<vec3<double>,
                vec3<double>,
                vec3<double>,
                vec3<double>,
                af::const_ref<double>,
                af::const_ref<vec3<double> > >((arg("initial_d1"),
                                                arg("initial_d2"),
                                                arg("initial_dn"),
                                                arg("initial_gp_offset"),
                                                arg("param_vals"),
                                                arg("param_axes"))))
      .def("d1", &PanelGroupCompose::d1)
      .def("d2", &PanelGroupCompose::d2)
      .def("origin", &PanelGroupCompose::origin)
      .def("derivatives_for_panel", &PanelGroupCompose::derivatives_for_panel);
  }

}}}  // namespace dials::refinement::boost_python
