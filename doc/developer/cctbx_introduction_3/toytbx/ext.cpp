#include <boost/python.hpp>
#include <cctype>

namespace toytbx { 
  namespace ext {

    static boost::python::list make_list(unsigned int n)
    {
      boost::python::list result;
      for(size_t i = 0; i < n; i++) {
	result.append(i);
      }
      return result;
    }

    void init_module()
    {
      using namespace boost::python;
      def("make_list", make_list, (arg("size")));
    }

  }
} // namespace toytbx::ext

BOOST_PYTHON_MODULE(toytbx_ext)
{
  toytbx::ext::init_module();
}
