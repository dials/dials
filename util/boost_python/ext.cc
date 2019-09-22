#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/util/scale_down_array.h>
#include <dials/util/masking.h>
#include <dials/util/export_mtz_helpers.h>

namespace dials { namespace util { namespace boost_python {

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_util_ext) {
    def("scale_down_array", &scale_down_array, (arg("image"), arg("scale_factor")));

    def("dials_u_to_mosflm", &dials_u_to_mosflm, (arg("dials_U"), arg("uc")));

    def("add_dials_batches",
        &add_dials_batches,
        (arg("mtz"),
         arg("dataset_id"),
         arg("image_range"),
         arg("batch_offset"),
         arg("wavelength"),
         arg("mosaic"),
         arg("phi_start"),
         arg("phi_range"),
         arg("cell_array"),
         arg("umat_array"),
         arg("panel_size"),
         arg("panel_distance"),
         arg("axis"),
         arg("s0n")));

    class_<ResolutionMaskGenerator>("ResolutionMaskGenerator", no_init)
      .def(init<const BeamBase &, const Panel &>())
      .def("apply", &ResolutionMaskGenerator::apply);
  }
}}}  // namespace dials::util::boost_python
