
#ifndef DIALS_ALGORITHMS_BOOST_PYTHON_RAY_PREDICTOR_WRAPPER_H
#define DIALS_ALGORITHMS_BOOST_PYTHON_RAY_PREDICTOR_WRAPPER_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spot_prediction/ray_predictor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename RayPredictorType>
  class_<RayPredictorType> ray_predictor_wrapper(const char *name) {
    
    // Useful typedefs
    typedef typename RayPredictorType::reflection_type reflection_type;
    typedef typename RayPredictorType::reflection_list_type reflection_list_type;

    // Typedef the different overloads for operator()
    reflection_list_type (RayPredictorType::*predict_single)(
      miller_index) const = &RayPredictorType::operator();
    reflection_list_type (RayPredictorType::*predict_array)(
      const flex_miller_index &) const = &RayPredictorType::operator();

    // Create and return the wrapper for the spot predictor object
    return class_ <RayPredictorType> (name, no_init)
      .def(init <vec3<double>,
                 vec3<double>,
                 mat3<double>,
                 vec2<double> > ((
        arg("s0"),
        arg("m2"),
        arg("UB"),
        arg("dphi"))))
      .def("__call__", predict_single, (
        arg("miller_index")))
      .def("__call__", predict_array, (
        arg("miller_indices")));
  }

}}} // namespace dials::algorithms::boost_python

#endif // DIALS_ALGORITHMS_BOOST_PYTHON_RAY_PREDICTOR_WRAPPER_H
