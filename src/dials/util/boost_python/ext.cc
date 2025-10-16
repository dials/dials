#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/util/scale_down_array.h>
#include <dials/util/masking.h>
#include <dials/util/export_mtz_helpers.h>
#include <dials/util/python_streambuf.h>

std::size_t dials::util::streambuf::default_buffer_size = 1024;
namespace dials { namespace util { namespace boost_python {
  struct python_streambuf_wrapper {
    typedef dials::util::streambuf wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("streambuf", no_init)
        .def(init<boost::python::object&, std::size_t>(
          (arg("python_file_obj"), arg("buffer_size") = 0)))
        .def_readwrite("default_buffer_size",
                       wt::default_buffer_size,
                       "The default size of the buffer sitting "
                       "between a Python file object and a C++ stream.");
    }
  };

  struct python_ostream_wrapper {
    typedef dials::util::ostream wt;

    static void wrap() {
      using namespace boost::python;
      class_<std::ostream, boost::noncopyable>("std_ostream", no_init);
      class_<wt, boost::noncopyable, bases<std::ostream> >("ostream", no_init)
        .def(init<boost::python::object&, std::size_t>(
          (arg("python_file_obj"), arg("buffer_size") = 0)));
    }
  };

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_util_ext) {
    def("scale_down_array", &scale_down_array, (arg("image"), arg("scale_factor")));

    def("ub_to_mosflm_u", &ub_to_mosflm_u, (arg("UB"), arg("uc")));

    class_<ResolutionMaskGenerator>("ResolutionMaskGenerator", no_init)
      .def(init<const BeamBase&, const Panel&>())
      .def("apply", &ResolutionMaskGenerator::apply);

    python_streambuf_wrapper::wrap();
    python_ostream_wrapper::wrap();
  }
}}}  // namespace dials::util::boost_python
