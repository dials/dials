
#ifndef DIALS_ALGORITHMS_BOOST_PYTHON_SPOT_PREDICTOR_WRAPPER_H
#define DIALS_ALGORITHMS_BOOST_PYTHON_SPOT_PREDICTOR_WRAPPER_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spot_prediction/spot_predictor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename SpotPredictorType>
  class_<SpotPredictorType> spot_predictor_wrapper(const char *name) {
    
    // Useful typedefs
    typedef typename SpotPredictorType::beam_type beam_type;
    typedef typename SpotPredictorType::detector_type detector_type;
    typedef typename SpotPredictorType::goniometer_type goniometer_type;
    typedef typename SpotPredictorType::scan_type scan_type;
    typedef typename SpotPredictorType::reflection_type reflection_type;
    typedef typename SpotPredictorType::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (SpotPredictorType::*predict_single)(
      miller_index) const = &SpotPredictorType::operator();
    reflection_list_type (SpotPredictorType::*predict_array)(
      const flex_miller_index &) const = &SpotPredictorType::operator();
    reflection_list_type (SpotPredictorType::*predict_generate)() = 
      &SpotPredictorType::operator();

    // Create and return the wrapper for the spot predictor object
    return class_ <SpotPredictorType> (name, no_init)
      .def(init <const beam_type&,
                 const detector_type&,
                 const goniometer_type&,
                 const scan_type&,
                 const cctbx::uctbx::unit_cell &,
                 const cctbx::sgtbx::space_group_type &,
                 mat3 <double>,
                 double> ((
        arg("beam"),
        arg("detector"),
        arg("goniometer"),
        arg("scan"),
        arg("unit_cell"),
        arg("space_group"),
        arg("ub_matrix"),
        arg("d_min"))))
      .def("__call__", predict_single, (
        arg("miller_index")))
      .def("__call__", predict_array, (
        arg("miller_indices")))
      .def("__call__", predict_generate);  
  }

}}} // namespace dials::algorithms::boost_python

#endif // DIALS_ALGORITHMS_BOOST_PYTHON_SPOT_PREDICTOR_WRAPPER_H