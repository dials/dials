
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform.h"

using namespace boost::python;

namespace dials { namespace integration { 
    
namespace boost_python {

void export_xds_transform() 
{
    void (XdsTransform::*calculate_single)(int, int, scitbx::af::tiny <int, 6>,
        scitbx::vec3 <double>, double) = &XdsTransform::calculate;

//    void (XdsTransform::*calculate_array)(const af::flex_tiny6_int &,
//        const af::flex_vec3_double &, const scitbx::af::flex_double &) = 
//            &XdsTransform::calculate;

    void (XdsTransform::*calculate_single_reflection)(int, Reflection &) = 
        &XdsTransform::calculate;

    void (XdsTransform::*calculate_array_reflection)(ReflectionList &) = 
        &XdsTransform::calculate;

    class_ <XdsTransform> ("XdsTransform")
        .def(init <XdsTransformGrid &,
                   const scitbx::af::flex_int &,
                   const scitbx::af::flex_int &,
                   const equipment::Detector&,
                   const equipment::Beam&,
                   const equipment::Goniometer&,
                   int > ((
                arg("xds_grid"),
                arg("image"),
                arg("mask"),
                arg("detector"),
                arg("beam"),
                arg("goniometer"),
                arg("n_div") = 5)))
        .def("calculate", calculate_single, (
                arg("reflection_index"),
                arg("mask_index"),
                arg("roi"),
                arg("s1"),
                arg("phi")))
//        .def("calculate", calculate_array, (
//                arg("roi"),
//                arg("s1"),
//                arg("phi")))
        .def("calculate", calculate_single_reflection, (
                arg("reflection_index"),
                arg("reflection")))
        .def("calculate", calculate_array_reflection, (
                arg("reflections")));
}

} // namespace = boost_python

}} // namespace = dials::integration
