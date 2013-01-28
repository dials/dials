
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../spot_predictor.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_spot_predictor()
{
    void (SpotPredictor::*predict_single)(cctbx::miller::index <>) =
        &SpotPredictor::predict;
    void (SpotPredictor::*predict_array)(const af::flex_miller_index &) =
        &SpotPredictor::predict;
    void (SpotPredictor::*predict_generate)() = &SpotPredictor::predict;
                
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
        .def("predict", predict_generate)        
        .add_property("miller_indices",
            &SpotPredictor::get_miller_indices)
        .add_property("rotation_angles",
            &SpotPredictor::get_rotation_angles)
        .add_property("beam_vectors",
            &SpotPredictor::get_beam_vectors)            
        .add_property("image_coordinates",
            &SpotPredictor::get_image_coordinates);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
