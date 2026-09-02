
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/tof/tof_mask_calculator.h>
#include <dials/algorithms/integration/tof/tof_integration.h>
#include <dials/algorithms/integration/tof/tof_profile_1d_ibix.h>
#include <dials/algorithms/integration/tof/tof_profile_3d_gutmann.h>
#include <dials/algorithms/scaling/tof/tof_scaling.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  std::shared_ptr<ProfileFitter> make_profile_fitter(
    object profile_1d_ibix_params_obj,
    object profile_3d_gutmann_params_obj) {
    // Only one profile fitting method allowed
    DIALS_ASSERT(!(!profile_1d_ibix_params_obj.is_none()
                   && !profile_3d_gutmann_params_obj.is_none()));

    if (!profile_1d_ibix_params_obj.is_none()) {
      return std::make_shared<Profile1DIBIXFitter>(
        extract<TOFProfile1DIBIXParams>(profile_1d_ibix_params_obj));
    }
    if (!profile_3d_gutmann_params_obj.is_none()) {
      return std::make_shared<Profile3DGutmannFitter>(
        extract<TOFProfile3DGutmannParams>(profile_3d_gutmann_params_obj));
    }
    return nullptr;
  }

  void integrate_reflection_table_wrapper(dials::af::reflection_table& reflection_table,
                                          dxtbx::model::Experiment& experiment,
                                          dxtbx::ImageSequence& data,
                                          object incident_params_obj,
                                          object absorption_params_obj,
                                          const bool& apply_lorentz,
                                          int n_threads,
                                          object profile_1d_ibix_params_obj,
                                          object profile_3d_gutmann_params_obj) {
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params;
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params;

    if (!incident_params_obj.is_none()) {
      incident_params =
        extract<dials_scaling::TOFIncidentSpectrumParams>(incident_params_obj);
    }
    if (!absorption_params_obj.is_none()) {
      absorption_params =
        extract<dials_scaling::TOFAbsorptionParams>(absorption_params_obj);
    }

    std::shared_ptr<ProfileFitter> profile_fitter =
      make_profile_fitter(profile_1d_ibix_params_obj, profile_3d_gutmann_params_obj);

    integrate_reflection_table(reflection_table,
                               experiment,
                               data,
                               incident_params,
                               absorption_params,
                               apply_lorentz,
                               n_threads,
                               profile_fitter);
  }

  void extract_correction_params(
    object incident_params_obj,
    object absorption_params_obj,
    boost::optional<dials_scaling::TOFIncidentSpectrumParams>& incident_params,
    boost::optional<dials_scaling::TOFAbsorptionParams>& absorption_params) {
    if (!incident_params_obj.is_none()) {
      incident_params =
        extract<dials_scaling::TOFIncidentSpectrumParams>(incident_params_obj);
    }
    if (!absorption_params_obj.is_none()) {
      absorption_params =
        extract<dials_scaling::TOFAbsorptionParams>(absorption_params_obj);
    }
  }

  boost::python::tuple calculate_line_profile_for_reflection_wrapper(
    dials::af::reflection_table& reflection,
    dxtbx::model::Experiment& experiment,
    dxtbx::ImageSequence& data,
    object incident_params_obj,
    object absorption_params_obj,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool& apply_lorentz_correction) {
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params;
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params;
    extract_correction_params(
      incident_params_obj, absorption_params_obj, incident_params, absorption_params);

    return calculate_line_profile_for_reflection(reflection,
                                                 experiment,
                                                 data,
                                                 incident_params,
                                                 absorption_params,
                                                 raw_projected_intensity_out,
                                                 projected_intensity_out,
                                                 projected_background_out,
                                                 tof_z_out,
                                                 apply_lorentz_correction);
  }

  boost::python::tuple calculate_line_profile_for_reflection_1d_wrapper(
    dials::af::reflection_table& reflection,
    dxtbx::model::Experiment& experiment,
    dxtbx::ImageSequence& data,
    object incident_params_obj,
    object absorption_params_obj,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    scitbx::af::shared<double> line_profile_out,
    const bool& apply_lorentz_correction,
    TOFProfile1DIBIXParams& profile_params_1d_ibix) {
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params;
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params;
    extract_correction_params(
      incident_params_obj, absorption_params_obj, incident_params, absorption_params);

    return calculate_line_profile_for_reflection(reflection,
                                                 experiment,
                                                 data,
                                                 incident_params,
                                                 absorption_params,
                                                 raw_projected_intensity_out,
                                                 projected_intensity_out,
                                                 projected_background_out,
                                                 tof_z_out,
                                                 line_profile_out,
                                                 apply_lorentz_correction,
                                                 profile_params_1d_ibix);
  }

  boost::python::tuple calculate_line_profile_for_reflection_3d_wrapper(
    dials::af::reflection_table& reflection,
    dxtbx::model::Experiment& experiment,
    dxtbx::ImageSequence& data,
    object incident_params_obj,
    object absorption_params_obj,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool& apply_lorentz_correction,
    TOFProfile3DGutmannParams& profile_params_3d_gutmann) {
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params;
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params;
    extract_correction_params(
      incident_params_obj, absorption_params_obj, incident_params, absorption_params);

    return calculate_line_profile_for_reflection_3d(reflection,
                                                    experiment,
                                                    data,
                                                    incident_params,
                                                    absorption_params,
                                                    raw_projected_intensity_out,
                                                    projected_intensity_out,
                                                    projected_background_out,
                                                    tof_z_out,
                                                    apply_lorentz_correction,
                                                    profile_params_3d_gutmann);
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
        &calculate_line_profile_for_reflection_wrapper,
        (arg("reflection"),
         arg("experiment"),
         arg("data"),
         arg("incident_params") = object(),
         arg("absorption_params") = object(),
         arg("raw_projected_intensity_out"),
         arg("projected_intensity_out"),
         arg("projected_background_out"),
         arg("tof_z_out"),
         arg("apply_lorentz_correction")));

    def("calculate_line_profile_for_reflection",
        &calculate_line_profile_for_reflection_1d_wrapper,
        (arg("reflection"),
         arg("experiment"),
         arg("data"),
         arg("incident_params"),
         arg("absorption_params"),
         arg("raw_projected_intensity_out"),
         arg("projected_intensity_out"),
         arg("projected_background_out"),
         arg("tof_z_out"),
         arg("line_profile_out"),
         arg("apply_lorentz_correction"),
         arg("profile_params_1d_ibix")));

    def("calculate_line_profile_for_reflection_3d",
        &calculate_line_profile_for_reflection_3d_wrapper,
        (arg("reflection"),
         arg("experiment"),
         arg("data"),
         arg("incident_params"),
         arg("absorption_params"),
         arg("raw_projected_intensity_out"),
         arg("projected_intensity_out"),
         arg("projected_background_out"),
         arg("tof_z_out"),
         arg("apply_lorentz_correction"),
         arg("profile_params_3d_gutmann")));
  }

}}}  // namespace dials::algorithms::boost_python
