
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../xds_transform.h"

using namespace boost::python;

namespace dials { namespace geometry { namespace algorithm { 
    
namespace boost_python {

void export_xds_transform() 
{
    class_ <XdsTransform> ("XdsTransform")
        .def(init <XdsTransformGrid,
                   scitbx::af::flex_int,
                   scitbx::vec3 <int>,
                   const equipment::Detector&,
                   const equipment::Beam&,
                   const equipment::Goniometer&,
                   scitbx::vec3 <int>,
                   int > ((
                arg("xds_grid"),
                arg("image"),
                arg("image_size"),
                arg("detector"),
                arg("beam"),
                arg("goniometer"),
                arg("roi_size"), 
                arg("n_div") = 5)))
        .def("calculate",
            &XdsTransform::calculate, (
                arg("reflection_index"),
                arg("xyz"),
                arg("s1"),
                arg("phi")));
}

} // namespace = boost_python

}}} // namespace = dials::geometry::algorithm
