
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/tof/tof_mask_calculator.h>
#include <dials/algorithms/integration/tof/tof_integration.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_algorithms_tof_integration_ext) {
    class_<TOFAbsorptionParams>("TOFAbsorptionParams", no_init)
      .def(init<double, double, double, double, double, double, double, double>());

    class_<TOFProfile1DParams>("TOFProfile1DParams", no_init)
      .def(init<double, double, double, double, double, double, double>());

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
        static_cast<void (*)(dials::af::reflection_table &,
                             dxtbx::model::Experiment &,
                             dxtbx::ImageSequence &,
                             bool)>(&integrate_reflection_table));

    def("integrate_reflection_table",
        static_cast<void (*)(dials::af::reflection_table &,
                             dxtbx::model::Experiment &,
                             dxtbx::ImageSequence &,
                             dxtbx::ImageSequence &,
                             dxtbx::ImageSequence &,
                             double,
                             double,
                             double,
                             bool)>(&integrate_reflection_table));

    def("integrate_reflection_table",
        static_cast<void (*)(dials::af::reflection_table &,
                             dxtbx::model::Experiment &,
                             dxtbx::ImageSequence &,
                             dxtbx::ImageSequence &,
                             dxtbx::ImageSequence &,
                             double,
                             double,
                             double,
                             TOFAbsorptionParams &,
                             bool)>(&integrate_reflection_table));

    def(
      "integrate_reflection_table_profile1d",
      static_cast<void (*)(dials::af::reflection_table &,
                           dxtbx::model::Experiment &,
                           dxtbx::ImageSequence &,
                           bool,
                           TOFProfile1DParams)>(&integrate_reflection_table_profile1d));
    def(
      "integrate_reflection_table_profile1d",
      static_cast<void (*)(dials::af::reflection_table &,
                           dxtbx::model::Experiment &,
                           dxtbx::ImageSequence &,
                           dxtbx::ImageSequence &,
                           dxtbx::ImageSequence &,
                           double,
                           double,
                           double,
                           bool,
                           TOFProfile1DParams)>(&integrate_reflection_table_profile1d));
    def("integrate_shoebox_profile1d",
        static_cast<boost::python::tuple (*)(dials::model::Shoebox<> &,
                                             dxtbx::model::Experiment &,
                                             dxtbx::ImageSequence &,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             scitbx::af::shared<double>,
                                             bool,
                                             TOFProfile1DParams)>(
          &integrate_shoebox_profile1d));
  }

}}}  // namespace dials::algorithms::boost_python
