#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/scaling/tof/tof_scaling.h>

namespace dials_scaling { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_tof_scaling_ext) {
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

    void (*extract_shoeboxes1)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes2)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               TOFIncidentSpectrumParams &,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes3)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               TOFIncidentSpectrumParams &,
                               TOFAbsorptionParams &,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;

    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes1);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes2);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes3);
  }

}}  // namespace dials_scaling::boost_python
