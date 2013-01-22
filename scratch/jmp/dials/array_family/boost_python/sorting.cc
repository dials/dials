
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../sorting.h"

using namespace boost::python;

namespace dials { namespace array_family { 

namespace boost_python {

void export_sorting() 
{
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_bool>);
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_int>);
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_long>);
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_size_t>);
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_float>);
    def("partial_sort_indices", 
        &partial_sort_indices_flex <scitbx::af::flex_double>);
}

} // namespace = boost_python

}} // namespace = dials::array_family
