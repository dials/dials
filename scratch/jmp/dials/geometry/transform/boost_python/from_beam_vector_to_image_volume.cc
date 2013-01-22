
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../from_beam_vector_to_image_volume.h"

using namespace boost::python;
using namespace dials::geometry::transform;

namespace dials { namespace geometry { namespace transform { 
    
namespace boost_python {

void export_from_beam_vector_to_image_volume() 
{
    // Apply to a single beam vector
    scitbx::vec3 <double> (FromBeamVectorToImageVolume::*apply_single)(
        scitbx::vec3 <double>, double) = &FromBeamVectorToImageVolume::apply;

    // Apply to array of beam vectors
    flex_vec3_double (FromBeamVectorToImageVolume::*apply_array)(
        const flex_vec3_double&, const scitbx::af::flex_double&, 
        scitbx::af::flex_bool&) = &FromBeamVectorToImageVolume::apply;

    class_ <FromBeamVectorToImageVolume> ("FromBeamVectorToImageVolume")
        .def(init <const equipment::Detector&,
                   const equipment::Goniometer&> ((
                arg("detector"), 
                arg("goniometer"))))
        .def("apply", apply_single, (
                arg("s1"),
                arg("phi")))
        .def("apply", apply_array, (
                arg("s1"),
                arg("phi"),
                arg("status")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::transform
