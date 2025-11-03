#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_function.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <dials/algorithms/symmetry/refstat/refstat.h>
#include <dials/algorithms/symmetry/refstat/merge_test.h>

namespace cctbx { namespace sgtbx { namespace refstat { namespace boost_python {

  namespace {
    struct nsp_wrapper {
      typedef named_space_group wt;

      static void wrap() {
        using namespace boost::python;
        return_value_policy<return_by_value> rbv;
        typedef return_internal_reference<> rir_t;
        class_<wt, bases<space_group>, std::auto_ptr<wt> >("named_space_group", no_init)
          .def(init<const space_group_symbols&>((arg("space_group_symbols"))))
          .add_property("name", make_function(&wt::get_name, rbv))
          .def("contains", &wt::contains_all);
        scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap(
          "shared_named_space_groups");
      }
    };

    template <typename FloatType>
    struct sae_wrapper {
      typedef extinction_element<FloatType> wt;

      static void wrap() {
        using namespace boost::python;
        typedef return_internal_reference<> rir_t;
        return_value_policy<return_by_value> rbv;
        class_<wt, std::auto_ptr<wt> >("extinction_element", no_init)
          .def_readonly("id", &wt::get_id)
          .def_readonly("name", &wt::name)
          .def_readonly("sumI", &wt::sumI)
          .def_readonly("sumS_sq", &wt::sumS_sq)
          .def_readonly("count", &wt::count)
          .def("rmx", &wt::get_rmx, rbv)
          .def("is_shadowed_by", &wt::is_shadowed_by)
          .def("shadowed_by_size", &wt::shadowed_by_cnt)
          .def("get_shadowed_by", &wt::get_shadowed_by, rir_t())
          .def("reset", &wt::reset)
          .def("__eq__", &wt::operator==);
        scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_extinctions");
      }
    };

    template <typename FloatType>
    struct sar_wrapper {
      typedef extinctions_registry<FloatType> wt;

      static void wrap() {
        using namespace boost::python;
        return_value_policy<return_by_value> rbv;
        typedef return_value_policy<copy_const_reference> ccr_t;
        typedef return_internal_reference<> rir_t;
        class_<wt>("extinctions_registry", no_init)
          .def(init<bool>((arg("init") = true)))
          .def("__len__", &wt::element_count)
          .def("__getitem__", &wt::get_element, rir_t())
          .def("__iter__", boost::python::iterator<wt, rir_t>())
          .def("get_extinctions", &wt::get_extinctions)
          .def("get_extinctions", &wt::get_extinctions_i)
          .def("find_sg", &wt::find_sg, rir_t())
          .def("sg_count", &wt::sg_count)
          .def("get_space_group", &wt::get_space_group, rir_t())
          .add_property("has_omp", &wt::has_openmp)
          .def("process",
               &wt::process,
               ((arg("indices"), arg("Is"), arg("sigs"), arg("scale") = 0)))
          .def("process_omp",
               &wt::process_omp,
               ((arg("indices"),
                 arg("Is"),
                 arg("sigs"),
                 arg("scale") = 0,
                 arg("threads_n") = -1)))
          .def("reset", &wt::reset)
          .def_readonly("ref_count", &wt::ref_count)
          .def_readonly("minI", &wt::minI)
          .def_readonly("maxI", &wt::maxI)
          .def_readonly("sumI", &wt::sumI)
          .def_readonly("sum_sig_sq", &wt::sum_sig_sq)
          .def_readonly("scale", &wt::scale);
      }
    };

    template <typename FloatType>
    struct merge_test_wrapper {
      struct wrap_merge_result {
        typedef typename merge_test_class<FloatType>::merge_stats_result wt;

        static void wrap() {
          using namespace boost::python;
          class_<wt, std::auto_ptr<wt> >("merge_stats_result", no_init)
            .def_readonly("inconsistent_count", &wt::inconsistent_count)
            .def_readonly("r_int", &wt::r_int);
        }
      };

      struct wrap_sysabs_test_result {
        typedef typename merge_test_class<FloatType>::sysabs_test_result wt;

        static void wrap() {
          using namespace boost::python;
          class_<wt, std::auto_ptr<wt> >("sysabs_test_result", no_init)
            .def_readonly("weak_count", &wt::weak_count)
            .def_readonly("weak_I_sum", &wt::weak_I_sum)
            .def_readonly("weak_sig_sq_sum", &wt::weak_sig_sq_sum)
            .def_readonly("strong_count", &wt::strong_count)
            .def_readonly("strong_I_sum", &wt::strong_I_sum)
            .def_readonly("strong_sig_sq_sum", &wt::strong_sig_sq_sum);
        }
      };

      struct wrap_test {
        typedef merge_test_class<FloatType> wt;

        static void wrap() {
          using namespace boost::python;
          class_<wt>("merge_test", no_init)
            .def(init<const af::shared<miller::index<> >&,
                      const af::shared<FloatType>&,
                      const af::shared<FloatType>&>(
              (arg("indices"), arg("Is"), arg("sigs"))))
            .def("merge_test", &wt::merge_test)
            .def("sysabs_test", &wt::sysabs_test);
        }
      };

      static void wrap() {
        wrap_merge_result::wrap();
        wrap_sysabs_test_result::wrap();
        wrap_test::wrap();
      }
    };

  }  // namespace

  void init_module() {
    nsp_wrapper::wrap();
    sae_wrapper<double>::wrap();
    sar_wrapper<double>::wrap();
    merge_test_wrapper<double>::wrap();
  }

}}}}  // namespace cctbx::sgtbx::refstat::boost_python

BOOST_PYTHON_MODULE(dials_algorithms_symmetry_refstat_ext) {
  cctbx::sgtbx::refstat::boost_python::init_module();
}
