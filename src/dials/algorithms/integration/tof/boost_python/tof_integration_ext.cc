
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/tof/tof_mask_calculator.h>
#include <dials/algorithms/integration/tof/tof_integration.h>
#include <dials/algorithms/integration/tof/tof_profile_1d_ibix.h>
#include <dials/algorithms/integration/tof/tof_profile_3d_gutmann.h>
#include <dials/algorithms/scaling/tof/tof_scaling.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void integrate_reflection_table_wrapper(dials::af::reflection_table& reflection_table,
                                          dxtbx::model::Experiment& experiment,
                                          dxtbx::ImageSequence& data,
                                          object incident_params_obj,
                                          object absorption_params_obj,
                                          const bool& apply_lorentz,
                                          int n_threads,
                                          object profile_1d_ibix_params_obj,
                                          object profile_3d_gutmann_params_obj) {
    boost::optional<TOFProfile1DIBIXParams> profile_1d_ibix_params;
    boost::optional<TOFProfile3DGutmannParams> profile_3d_gutmann_params;

    if (!profile_1d_ibix_params_obj.is_none()) {
      profile_1d_ibix_params =
        extract<TOFProfile1DIBIXParams>(profile_1d_ibix_params_obj);
    }

    if (!profile_3d_gutmann_params_obj.is_none()) {
      profile_3d_gutmann_params =
        extract<TOFProfile3DGutmannParams>(profile_3d_gutmann_params_obj);
    }

    if (absorption_params_obj.is_none() && incident_params_obj.is_none()) {
      integrate_reflection_table(reflection_table,
                                 experiment,
                                 data,
                                 apply_lorentz,
                                 n_threads,
                                 profile_1d_ibix_params,
                                 profile_3d_gutmann_params);

      return;
    }

    if (!incident_params_obj.is_none()) {
      dials_scaling::TOFIncidentSpectrumParams incident_params =
        extract<dials_scaling::TOFIncidentSpectrumParams>(incident_params_obj);

      if (!absorption_params_obj.is_none()) {
        dials_scaling::TOFAbsorptionParams absorption_params =
          extract<dials_scaling::TOFAbsorptionParams>(absorption_params_obj);

        integrate_reflection_table(reflection_table,
                                   experiment,
                                   data,
                                   incident_params,
                                   absorption_params,
                                   apply_lorentz,
                                   n_threads,
                                   profile_1d_ibix_params,
                                   profile_3d_gutmann_params);
      }

      else {
        integrate_reflection_table(reflection_table,
                                   experiment,
                                   data,
                                   incident_params,
                                   apply_lorentz,
                                   n_threads,
                                   profile_1d_ibix_params,
                                   profile_3d_gutmann_params);
      }
    }
  }

  BOOST_PYTHON_MODULE(dials_algorithms_tof_integration_ext) {
    class_<TOFProfile1DIBIXParams>("TOFProfile1DIBIXParams", no_init)
      .def(
        init<double, double, double, double, double, double, double, int, bool, bool>())
      .def_readwrite("A", &TOFProfile1DIBIXParams::A)
      .def_readwrite("alpha", &TOFProfile1DIBIXParams::alpha)
      .def_readwrite("alpha_min", &TOFProfile1DIBIXParams::alpha_min)
      .def_readwrite("alpha_max", &TOFProfile1DIBIXParams::alpha_max)
      .def_readwrite("beta", &TOFProfile1DIBIXParams::beta)
      .def_readwrite("beta_min", &TOFProfile1DIBIXParams::beta_min)
      .def_readwrite("beta_max", &TOFProfile1DIBIXParams::beta_max)
      .def_readwrite("n_restarts", &TOFProfile1DIBIXParams::n_restarts)
      .def_readwrite("optimize_profile", &TOFProfile1DIBIXParams::optimize_profile)
      .def_readwrite("show_profile_failures",
                     &TOFProfile1DIBIXParams::show_profile_failures);

    class_<TOFProfile3DGutmannParams>("TOFProfile3DGutmannParams", no_init)
      .def(
        init<double, double, double, double, double, double, int, bool, bool, bool>())
      .def_readwrite("alpha", &TOFProfile3DGutmannParams::alpha)
      .def_readwrite("alpha_min", &TOFProfile3DGutmannParams::alpha_min)
      .def_readwrite("alpha_max", &TOFProfile3DGutmannParams::alpha_max)
      .def_readwrite("beta", &TOFProfile3DGutmannParams::beta)
      .def_readwrite("beta_min", &TOFProfile3DGutmannParams::beta_min)
      .def_readwrite("beta_max", &TOFProfile3DGutmannParams::beta_max)
      .def_readwrite("n_restarts", &TOFProfile3DGutmannParams::n_restarts)
      .def_readwrite("optimize_profile", &TOFProfile3DGutmannParams::optimize_profile)
      .def_readwrite("use_central_diff", &TOFProfile3DGutmannParams::use_central_diff)
      .def_readwrite("show_profile_failures",
                     &TOFProfile3DGutmannParams::show_profile_failures);

    def("tof_calculate_ellipse_shoebox_mask",
        &tof_calculate_ellipse_shoebox_mask,
        (arg("reflection_table"),
         arg("experiment"),
         arg("n_threads") = 1,
         arg("scale") = 1));

    def("tof_calculate_seed_skewness_shoebox_mask",
        &tof_calculate_seed_skewness_shoebox_mask,
        (arg("reflection_table"),
         arg("experiment"),
         arg("d_skewness_threshold"),
         arg("min_iterations"),
         arg("n_threads") = 1));

    def("integrate_reflection_table",
        &integrate_reflection_table_wrapper,
        (arg("reflection_table"),
         arg("experiment"),
         arg("data"),
         arg("incident_params"),
         arg("absorption_params"),
         arg("apply_lorentz_correction"),
         arg("n_threads"),
         arg("profile_1d_ibix_params") = object()));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table&,
                                             dxtbx::model::Experiment&,
                                             dxtbx::ImageSequence&,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool&)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection_3d",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table&,
                                             dxtbx::model::Experiment&,
                                             dxtbx::ImageSequence&,
                                             scitbx::af::shared<vec3<double> >,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool&,
                                             TOFProfile3DGutmannParams&)>(
          &calculate_line_profile_for_reflection_3d));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table&,
                                             dxtbx::model::Experiment&,
                                             dxtbx::ImageSequence&,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool&,
                                             TOFProfile1DIBIXParams&)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(
          dials::af::reflection_table&,
          dxtbx::model::Experiment&,
          dxtbx::ImageSequence&,
          const dials_scaling::TOFIncidentSpectrumParams&,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          const bool&)>(&calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(
          dials::af::reflection_table&,
          dxtbx::model::Experiment&,
          dxtbx::ImageSequence&,
          const dials_scaling::TOFIncidentSpectrumParams&,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          const bool&,
          TOFProfile1DIBIXParams&)>(&calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(
          dials::af::reflection_table&,
          dxtbx::model::Experiment&,
          dxtbx::ImageSequence&,
          const dials_scaling::TOFIncidentSpectrumParams&,
          const dials_scaling::TOFAbsorptionParams&,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          const bool&)>(&calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(
          dials::af::reflection_table&,
          dxtbx::model::Experiment&,
          dxtbx::ImageSequence&,
          const dials_scaling::TOFIncidentSpectrumParams&,
          const dials_scaling::TOFAbsorptionParams&,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          scitbx::af::shared<double>,
          const bool&,
          TOFProfile1DIBIXParams&)>(&calculate_line_profile_for_reflection));
  }

}}}  // namespace dials::algorithms::boost_python
