#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <cctype>

namespace x2tbx { 
  namespace ext {

    // make a python list

    static boost::python::list make_list(size_t n)
    {
      boost::python::list result;
      for(size_t i = 0; i < n; i++) {
	result.append(i);
      }
      return result;
    }

    // make a flex array (much more flexible)

    static scitbx::af::shared<int> make_flex(size_t n)
    {
      scitbx::af::shared<int> result;
      for(size_t i = 0; i < n; i++) {
	result.push_back(i);
      }
      return result;
    }

    // using flex arrays

    static int sum(scitbx::af::shared<int> array)
    {
      int result = 0;
      for (size_t i = 0; i < array.size(); i++) {
	result += array[i];
      }
      return result;
    }

    void init_module()
    {
      using namespace boost::python;
      def("make_list", make_list, (arg("size")));
      def("make_flex", make_flex, (arg("size")));
      def("sum", sum, (arg("array")));
    }

  }
} // namespace x2tbx::ext

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::ext::init_module();
}
