#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/scaling/tof/tof_scaling.h>

namespace dials_scaling { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_tof_scaling_ext) {
    class_<TOFCorrectionsData>("TOFCorrectionsData", no_init)
      .def(init<double,
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double>());

    void (*extract_shoeboxes1)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes2)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               dxtbx::ImageSequence &,
                               dxtbx::ImageSequence &,
                               TOFCorrectionsData &,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;
    void (*extract_shoeboxes3)(dials::af::reflection_table &,
                               dxtbx::model::Experiment &,
                               dxtbx::ImageSequence &,
                               dxtbx::ImageSequence &,
                               dxtbx::ImageSequence &,
                               double,
                               double,
                               double,
                               bool) = &tof_extract_shoeboxes_to_reflection_table;

    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes1);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes2);
    def("tof_extract_shoeboxes_to_reflection_table", extract_shoeboxes3);
  }

}}  // namespace dials_scaling::boost_python
