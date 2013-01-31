
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/format.hpp>
#include <string>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include "../reflection.h"

using namespace boost::python;

namespace dials { 

namespace boost_python {

std::string reflection_to_string(const Reflection &reflection) {
    boost::format fmt(
        "Reflection:\n"
        "    miller index:   (%1%, %2%, %3%)\n"
        "    rotation angle: %4%\n"
        "    beam vector:    (%5%, %6%, %7%)\n"
        "    image coord:    (%8%, %9%, %10%)\n"
        "    ROI:            (%11%, %12%, %13%, %14%, %15%, %16%)\n"
        "    Mask index:     %17%\n"
        "    Background:     %18%");
        
    fmt % reflection.get_miller_index()[0];
    fmt % reflection.get_miller_index()[1];
    fmt % reflection.get_miller_index()[2];
    fmt % reflection.get_rotation_angle();
    fmt % reflection.get_beam_vector()[0];
    fmt % reflection.get_beam_vector()[1];
    fmt % reflection.get_beam_vector()[2];
    fmt % reflection.get_image_coord()[0];
    fmt % reflection.get_image_coord()[1];
    fmt % reflection.get_image_coord()[2];
    fmt % reflection.get_region_of_interest()[0];
    fmt % reflection.get_region_of_interest()[1];
    fmt % reflection.get_region_of_interest()[2];
    fmt % reflection.get_region_of_interest()[3];
    fmt % reflection.get_region_of_interest()[4];
    fmt % reflection.get_region_of_interest()[5];
    fmt % reflection.get_mask_index();
    fmt % reflection.get_background_intensity();
    return fmt.str();
}

void export_reflection()
{
    class_<Reflection> ("Reflection")
        .def(init <cctbx::miller::index <>,
                   double,
                   scitbx::vec3 <double>,
                   scitbx::vec3 <double> > ((
                        arg("miller_index"),
                        arg("rotation_angle"),
                        arg("beam_vector"),
                        arg("image_coord"))))
        .add_property("miller_index", 
            &Reflection::get_miller_index,
            &Reflection::set_miller_index)
        .add_property("rotation_angle", 
            &Reflection::get_rotation_angle,
            &Reflection::set_rotation_angle)
        .add_property("beam_vector", 
            &Reflection::get_beam_vector,
            &Reflection::set_beam_vector)
        .add_property("image_coord",
            &Reflection::get_image_coord,
            &Reflection::set_image_coord)
        .add_property("region_of_interest",
            &Reflection::get_region_of_interest,
            &Reflection::set_region_of_interest)
        .add_property("mask_index",
            &Reflection::get_mask_index,
            &Reflection::set_mask_index)
        .add_property("background_intensity",
            &Reflection::get_background_intensity,
            &Reflection::set_background_intensity)
        .def("is_zero", &Reflection::is_zero)
        .def("__str__", &reflection_to_string);          
        
    scitbx::af::boost_python::flex_wrapper <Reflection>::plain("ReflectionList");        
    //class_<ReflectionList>("ReflectionList")
    //    .def(map_indexing_suite<ReflectionList>());        
}

} // namespace = boost_python

} // namespace = dials
