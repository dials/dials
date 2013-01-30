
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../spot_predictor.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_spot_predictor()
{
    scitbx::vec2 <Reflection> (SpotPredictor::*predict_single)(
        cctbx::miller::index <>) const = &SpotPredictor::predict;
    scitbx::af::shared <Reflection> (SpotPredictor::*predict_array)(
        const af::flex_miller_index &) const = &SpotPredictor::predict;
    scitbx::af::shared <Reflection> (SpotPredictor::*predict_generate)() = 
        &SpotPredictor::predict;
                
    class_ <SpotPredictor> ("SpotPredictor", no_init)
        .def(init <const equipment::Beam &,
                   const equipment::Detector &,
                   const equipment::Goniometer &,
                   const cctbx::uctbx::unit_cell &,
                   const cctbx::sgtbx::space_group_type &,
                   scitbx::mat3 <double>,
                   double> ((
            arg("beam"),
            arg("detector"),
            arg("goniometer"),
            arg("unit_cell"),
            arg("space_group_type"),
            arg("ub_matrix"),
            arg("d_min"))))
        .def("predict", predict_single, (
            arg("miller_index")))
        .def("predict", predict_array, (
            arg("miller_indices")))
        .def("predict", predict_generate);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
