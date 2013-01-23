#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace dials { namespace af {

namespace boost_python {

void export_flex_tiny()
{
    using namespace boost::python;
    scitbx::af::boost_python::flex_wrapper <scitbx::af::tiny <int, 6> >
        ::plain("tiny6_int");
    scitbx::af::boost_python::flex_wrapper <scitbx::af::tiny <double, 6> >
        ::plain("tiny6_double");
}

} // namespace boost_python

}} // namespace dials::af
