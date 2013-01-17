
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../index_generator.h"

using namespace boost::python;

namespace dials { namespace spot_prediction { 

namespace boost_python {

void export_index_generator()
{
    class_ <IndexGenerator> ("IndexGenerator")
        .def(init <cctbx::uctbx::unit_cell const&,
                   cctbx::sgtbx::space_group_type const&,
                   bool,
                   double> ((
            arg("unit_cell"),
            arg("space_group_type"),
            arg("anomalous_flag"),
            arg("resolution_d_min"))))
        .def("next", &IndexGenerator::next)
        .def("to_array", &IndexGenerator::to_array);
}

} // namespace = boost_python

}} // namespace = dials::spot_prediction
