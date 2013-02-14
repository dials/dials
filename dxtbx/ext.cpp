#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <boost_adaptbx/python_streambuf.h>
#include <cctype>
#include <fstream>

namespace dxtbx {
  namespace ext {

    scitbx::af::shared<int>
    read_uint16(boost_adaptbx::python::streambuf & input,
		size_t count)
    {
      scitbx::af::shared<int> result;
      boost_adaptbx::python::streambuf::istream is(input);
      unsigned short * data = new unsigned short[count];

      is.read((char *) data, count * sizeof(unsigned short));

      for (size_t j = 0; j < count; j++) {
	result.push_back((int) data[j]);
      }

      delete[](data);

      return result;
    }

    void init_module()
    {
      using namespace boost::python;
      def("read_uint16", read_uint16, (arg("file"), arg("count")));
    }
  }
} //namespace dxtbx::ext

BOOST_PYTHON_MODULE(dxtbx_ext)
{
  dxtbx::ext::init_module();
}
