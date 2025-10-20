
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/tof/tof_mask_calculator.h>
#include <dials/algorithms/integration/tof/tof_integration.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void integrate_reflection_table_wrapper(dials::af::reflection_table &reflection_table,
                                          dxtbx::model::Experiment &experiment,
                                          dxtbx::ImageSequence &data,
                                          object incident_params_obj,
                                          object absorption_params_obj,
                                          const bool &apply_lorentz,
                                          int n_threads,
                                          object profile1d_params_obj,
                                          object profile3d_params_obj) {
    boost::optional<TOFProfile1DParams> profile1d_params;
    boost::optional<TOFProfile3DParams> profile3d_params;

    if (!profile1d_params_obj.is_none()) {
      profile1d_params = extract<TOFProfile1DParams>(profile1d_params_obj);
    }

    if (!profile3d_params_obj.is_none()) {
      profile3d_params = extract<TOFProfile3DParams>(profile3d_params_obj);
    }

    if (absorption_params_obj.is_none() && incident_params_obj.is_none()) {
      integrate_reflection_table(reflection_table,
                                 experiment,
                                 data,
                                 apply_lorentz,
                                 n_threads,
                                 profile1d_params,
                                 profile3d_params);

      return;
    }

    if (!incident_params_obj.is_none()) {
      TOFIncidentSpectrumParams incident_params =
        extract<TOFIncidentSpectrumParams>(incident_params_obj);

      if (!absorption_params_obj.is_none()) {
        TOFAbsorptionParams absorption_params =
          extract<TOFAbsorptionParams>(absorption_params_obj);

        integrate_reflection_table(reflection_table,
                                   experiment,
                                   data,
                                   incident_params,
                                   absorption_params,
                                   apply_lorentz,
                                   n_threads,
                                   profile1d_params,
                                   profile3d_params);
      }

      else {
        integrate_reflection_table(reflection_table,
                                   experiment,
                                   data,
                                   incident_params,
                                   apply_lorentz,
                                   n_threads,
                                   profile1d_params,
                                   profile3d_params);
      }
    }
  }
  void test(scitbx::af::versa<double, scitbx::af::c_grid<3> > ye) {}

  BOOST_PYTHON_MODULE(dials_algorithms_tof_integration_ext) {
    class_<TOFAbsorptionParams>("TOFAbsorptionParams", no_init)
      .def(init<double, double, double, double, double, double, double, double>());

    class_<TOFIncidentSpectrumParams>("TOFIncidentSpectrumParams", no_init)
      .def(init<std::shared_ptr<ImageSequence>,
                std::shared_ptr<ImageSequence>,
                double,
                double,
                double>())
      .def_readwrite("incident_data", &TOFIncidentSpectrumParams::incident_data)
      .def_readwrite("empty_data", &TOFIncidentSpectrumParams::empty_data)
      .def_readwrite("sample_proton_charge",
                     &TOFIncidentSpectrumParams::sample_proton_charge)
      .def_readwrite("incident_proton_charge",
                     &TOFIncidentSpectrumParams::incident_proton_charge)
      .def_readwrite("empty_proton_charge",
                     &TOFIncidentSpectrumParams::empty_proton_charge);

    class_<TOFProfile1DParams>("TOFProfile1DParams", no_init)
      .def(init<double, double, double, double, double, double, double, int>());

    class_<TOFProfile3DParams>("TOFProfile3DParams", no_init)
      .def(init<double, double, double, double, double, double, int>());

    def("tof_calculate_ellipse_shoebox_mask",
        &tof_calculate_ellipse_shoebox_mask,
        (arg("reflection_table"), arg("experiment")));

    def("tof_calculate_seed_skewness_shoebox_mask",
        &tof_calculate_seed_skewness_shoebox_mask,
        (arg("reflection_table"),
         arg("experiment"),
         arg("d_skewness_threshold"),
         arg("min_iterations")));

    def("extract_shoeboxes_to_reflection_table",
        &extract_shoeboxes_to_reflection_table,
        (arg("reflection_table"), arg("experiment"), arg("data")));

    def("integrate_reflection_table",
        &integrate_reflection_table_wrapper,
        (arg("reflection_table"),
         arg("experiment"),
         arg("data"),
         arg("incident_params"),
         arg("absorption_params"),
         arg("apply_lorentz_correction"),
         arg("n_threads"),
         arg("profile1d_params") = object()));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection_3d",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &,
                                             const TOFProfile3DParams &)>(
          &calculate_line_profile_for_reflection_3d));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &,
                                             const TOFProfile1DParams &)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             const TOFIncidentSpectrumParams &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             const TOFIncidentSpectrumParams &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &,
                                             const TOFProfile1DParams &)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             const TOFIncidentSpectrumParams &,
                                             const TOFAbsorptionParams &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &)>(
          &calculate_line_profile_for_reflection));

    def("calculate_line_profile_for_reflection",
        static_cast<boost::python::tuple (*)(dials::af::reflection_table &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             const TOFIncidentSpectrumParams &,
                                             const TOFAbsorptionParams &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             const bool &,
                                             const TOFProfile1DParams &)>(
          &calculate_line_profile_for_reflection));
  }

}}}  // namespace dials::algorithms::boost_python
