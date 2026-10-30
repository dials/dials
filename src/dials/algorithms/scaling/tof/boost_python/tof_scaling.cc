#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/scaling/tof/tof_scaling.h>
#include <dials/algorithms/scaling/tof/lambda_norm.h>

namespace dials_scaling { namespace boost_python {

  using namespace boost::python;

  void export_lambda_norm() {
    class_<LambdaNormVarianceModel>("LambdaNormVarianceModel")
      .def("fit", &LambdaNormVarianceModel::fit)
      .def("b", &LambdaNormVarianceModel::b)
      .def("variance", &LambdaNormVarianceModel::variance)
      .def("lambda_edges", &LambdaNormVarianceModel::lambda_edges)
      .def("b_values", &LambdaNormVarianceModel::b_values)
      .def("counting_share", &LambdaNormVarianceModel::counting_share);

    class_<LambdaNormSamples>("LambdaNormSamples", no_init)
      .add_property(
        "index",
        make_getter(&LambdaNormSamples::index, return_value_policy<return_by_value>()))
      .add_property(
        "rho",
        make_getter(&LambdaNormSamples::rho, return_value_policy<return_by_value>()))
      .add_property(
        "u",
        make_getter(&LambdaNormSamples::u, return_value_policy<return_by_value>()));

    def("leave_one_out_samples", &leave_one_out_samples);

    class_<GPRWavelengthScale>("GPRWavelengthScale", no_init)
      .def(init<std::size_t, double, bool, double, double>((arg("nbins"),
                                                            arg("lengthscale"),
                                                            arg("matern32"),
                                                            arg("lambda_min"),
                                                            arg("lambda_max"))))
      .def("fit", &GPRWavelengthScale::fit)
      .def("scale", &GPRWavelengthScale::scale)
      .def("uncertainty", &GPRWavelengthScale::uncertainty)
      .def("active", &GPRWavelengthScale::active)
      .def("message", &GPRWavelengthScale::message)
      .def("lambda_min", &GPRWavelengthScale::lambda_min)
      .def("lambda_max", &GPRWavelengthScale::lambda_max)
      .def("lambda_ref", &GPRWavelengthScale::lambda_ref)
      .def("lengthscale", &GPRWavelengthScale::lengthscale)
      .def("sigma_f", &GPRWavelengthScale::sigma_f)
      .def("noise_inflation", &GPRWavelengthScale::noise_inflation)
      .def("mean_noise", &GPRWavelengthScale::mean_noise)
      .def("log_marginal_likelihood", &GPRWavelengthScale::log_marginal_likelihood)
      .def("n_samples", &GPRWavelengthScale::n_samples)
      .def("lengthscale_at_lower_limit",
           &GPRWavelengthScale::lengthscale_at_lower_limit)
      .def("lengthscale_at_upper_limit",
           &GPRWavelengthScale::lengthscale_at_upper_limit)
      .def("noise_inflation_at_upper_limit",
           &GPRWavelengthScale::noise_inflation_at_upper_limit)
      .def("training_lambda", &GPRWavelengthScale::training_lambda)
      .def("training_scale", &GPRWavelengthScale::training_scale)
      .def("training_sigma", &GPRWavelengthScale::training_sigma);
  }

  BOOST_PYTHON_MODULE(dials_tof_scaling_ext) {
    export_lambda_norm();

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

    void (*extract_shoeboxes1)(dials::af::reflection_table&,
                               dxtbx::model::Experiment&,
                               dxtbx::ImageSequence&,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes2)(dials::af::reflection_table&,
                               dxtbx::model::Experiment&,
                               dxtbx::ImageSequence&,
                               TOFIncidentSpectrumParams&,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes3)(dials::af::reflection_table&,
                               dxtbx::model::Experiment&,
                               dxtbx::ImageSequence&,
                               TOFIncidentSpectrumParams&,
                               TOFAbsorptionParams&,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;

    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes1);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes2);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes3);
  }

}}  // namespace dials_scaling::boost_python
