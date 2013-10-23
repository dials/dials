#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/serialization/single_buffered.h>
#include <scitbx/matrix/transpose_multiply.h>
#include <scitbx/math/utils.h>
#include <scitbx/array_family/tiny_types.h>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/format.hpp>
#include <scitbx/array_family/boost_python/flex_helpers.h>

namespace scitbx { namespace serialization { namespace single_buffered {

  inline
  char* to_string(char* start, scitbx::af::int2 const& value)
  {
    return
      to_string(to_string(start, value[0]), value[1]);
  }

  template <>
  struct from_string<scitbx::af::int2>
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<2;i++) {
        from_string<int> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    scitbx::af::int2 value;
    const char* end;
  };

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  template <>
  struct flex_default_element<scitbx::af::int2>
  {
    static scitbx::af::int2
    get() { return scitbx::af::int2(0,0); }
  };

}}}

namespace dials { namespace af {
namespace {

  scitbx::af::flex<scitbx::af::int2>::type*
  join(
    scitbx::af::const_ref<int> const& a,
    scitbx::af::const_ref<int> const& b)
  {
    SCITBX_ASSERT(a.size() == b.size());
    scitbx::af::shared<scitbx::af::int2> result((scitbx::af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(scitbx::af::int2(a[i],b[i]));
    }
    return new scitbx::af::flex<scitbx::af::int2>::type(result, result.size());
  }

  scitbx::af::flex<scitbx::af::int2>::type*
  from_int(
    scitbx::af::const_ref<int> const& x)
  {
    SCITBX_ASSERT(x.size() % 2 == 0);
    std::size_t result_size = x.size() / 2;
    scitbx::af::shared<scitbx::af::int2> result(result_size);
    const int* d = x.begin();
    for(std::size_t i=0;i<result_size;i++) {
      for (std::size_t j=0;j<2;++j) {
        result[i][j] = *d++;
      }
    }
    return new scitbx::af::flex<scitbx::af::int2>::type(result, result.size());
  }

  scitbx::af::flex_int
  as_int(scitbx::af::flex<scitbx::af::int2>::type const& a)
  {
    SCITBX_ASSERT(a.accessor().is_trivial_1d());
    scitbx::af::flex_int result(a.size()*2, scitbx::af::init_functor_null<int>());
    int* r = result.begin();
    scitbx::af::const_ref<scitbx::af::int2> a_ref = a.const_ref().as_1d();
    for(std::size_t i=0;i<a_ref.size();i++) {
      for(std::size_t j=0;j<2;j++) {
        *r++ = a_ref[i][j];
      }
    }
    return result;
  }

} // namespace <anonymous>

namespace boost_python {



  void export_flex_int2()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef scitbx::af::boost_python::flex_wrapper<scitbx::af::int2> f_w;
    f_w::plain("int2")
      .def_pickle(scitbx::af::boost_python::flex_pickle_single_buffered<
        scitbx::af::int2,
        2*scitbx::af::boost_python::pickle_size_per_element<int>::value>())
      .def("__init__", make_constructor(join))
      .def("__init__", make_constructor(from_int))
      .def("as_int", as_int);
    ;
  }

}}} // namespace scitbx::af::boost_python
