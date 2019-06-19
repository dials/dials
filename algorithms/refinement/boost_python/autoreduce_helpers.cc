#include <dials/array_family/flex_table.h>
#include <dials/array_family/reflection_table.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <algorithm>
#include <numeric>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <iterator>

using namespace boost::python;
using namespace dials::af;

// Helper classes to speed up calculations in autoreduce.py

namespace dials { namespace refinement { namespace boost_python {
  /*
    Helper class used by the AutoReduce._panel_gp_surplus_reflections function in
    autoreduce.py
  */
  struct pg_surpl_iter {
    template <typename IntType>
    pg_surpl_iter(const scitbx::af::shared<IntType>& ref_ids,
                  const scitbx::af::shared<std::size_t>& panel_id,
                  const boost::python::list pnl_ids,
                  const scitbx::af::shared<std::size_t> exp_ids,
                  const int cutoff) {
      // Perform python list conversion once; improves read access
      scitbx::af::shared<std::size_t> pnl_ids_flex;
      for (int ii = 0; ii < len(pnl_ids); ++ii) {
        pnl_ids_flex.push_back(boost::python::extract<std::size_t>(pnl_ids[ii]));
      }

      std::size_t count = 0;

      // Work with raw pointers instead of overloadable objects
      const std::size_t* ptr_exp_ids = exp_ids.begin();
      const std::size_t* ptr_panel_id = panel_id.begin();

      // Predefine loop variables
      std::size_t exp_it, pnl;

      std::pair<const IntType*, const IntType*> refIDRange;
      std::pair<const std::size_t*, const std::size_t*> panelIDRange;
      unsigned int data_offset = 0;
      unsigned int data_range = 0;

      /* Search for the beginning and end indices for the index values and track their
         pointers. Use these pointers to narow the next search margin and count the
         range of elements that match. These values indicate the count and can be added
         to a running total.
      */
      for (exp_it = 0; exp_it < exp_ids.size();
           ++exp_it) {  // Iterate through experiment IDs
        // Get the pointer indices to the beginning and end of the values with the
        // requested experimental ID. Performs search in ~2logN
        refIDRange =
          std::equal_range(ref_ids.begin(), ref_ids.end(), ptr_exp_ids[exp_it]);

        // Calculate pointer offset for start and end of the data to examine.
        data_range = refIDRange.second - refIDRange.first;
        data_offset = refIDRange.first - ref_ids.begin();

        // Taking the previous pointer positions, and perfoming a search within this
        // range for the matching panel IDs, we then count the pointer range as all
        // valid entries and add to the accumlator.
        for (pnl = 0; pnl < pnl_ids_flex.size(); ++pnl) {  // Iterate through panel IDs
          panelIDRange = std::equal_range(ptr_panel_id + data_offset,
                                          ptr_panel_id + data_offset + data_range,
                                          pnl_ids_flex[pnl]);
          count += (unsigned int)(panelIDRange.second - panelIDRange.first);
        }
      }
      result = count - cutoff;
    }
    int result;
  };

  /*
    Helper class used by the AutoReduce._unit_cell_surplus_reflections function in
    autoreduce.py
  */
  struct uc_surpl_iter {
    template <typename IntType>
    uc_surpl_iter(const scitbx::af::shared<IntType>& ref_ids,
                  const scitbx::af::shared<cctbx::miller::index<> >& ref_miller,
                  const scitbx::af::shared<std::size_t>& exp_ids,
                  const scitbx::af::shared<mat3<double> > F_dbdp) {
      // Expose pointers for experiment flex arrays
      const std::size_t* ptr_exp_ids = exp_ids.begin();

      std::size_t exp_it;  // ref_id;
      std::vector<std::pair<const IntType*, const IntType*> > refIDRange(
        exp_ids.size());

      std::vector<unsigned int> data_range(exp_ids.size());
      std::vector<unsigned int> data_offset(exp_ids.size());

      // Perform subselection of data ranges
      for (exp_it = 0; exp_it < exp_ids.size(); ++exp_it) {
        refIDRange[exp_it] =
          std::equal_range(ref_ids.begin(), ref_ids.end(), ptr_exp_ids[exp_it]);
        data_range[exp_it] = refIDRange[exp_it].second - refIDRange[exp_it].first;
        data_offset[exp_it] = refIDRange[exp_it].first - ref_ids.begin();
      }

      // Pre-initialise loop variables and array for storing norms
      int ii, jj;
      scitbx::af::shared<double> mul_norm;

      // Calculate L-2 norm of Matrix*Vec operation for each matrix in F_dbdp and each
      // vec3 miller index across all exp_ids
      for (jj = 0; jj < F_dbdp.size(); ++jj) {
        for (unsigned int exp_id_range = 0; exp_id_range < data_range.size();
             ++exp_id_range) {
          for (ii = data_offset[exp_id_range];
               ii < data_offset[exp_id_range] + data_range[exp_id_range];
               ++ii) {
            vec3<double> v3(ref_miller[ii]);
            mul_norm.push_back((F_dbdp[jj] * v3).length());
          }
        }
      }

      // Sum the total offsets to determine each range for splitting of mul_norm to
      // nref param values with F_dbdp.size() number of partitions
      int sum_offsets = std::accumulate(data_range.begin(), data_range.end(), 0);
      scitbx::af::shared<int> nref_each_param;
      for (ii = 0; ii < F_dbdp.size(); ++ii) {
        // Count the number of elements greater than 0 from the L-2 norms.
        nref_each_param.push_back(
          std::count_if(mul_norm.begin() + ii * sum_offsets,
                        mul_norm.begin() + (ii + 1) * sum_offsets,
                        gtZero));
      }

      // Return the value of the smallest element
      scitbx::af::shared<int>::iterator it =
        std::min_element(nref_each_param.begin(), nref_each_param.end());
      result = *(nref_each_param.begin() + std::distance(nref_each_param.begin(), it));
    }
    int result;
    static bool gtZero(double x) {
      return ((x > 0.0) == 1);
    }  // Used for comparison operation
  };

  struct surpl_iter {
    template <typename IntType>
    surpl_iter(scitbx::af::shared<IntType> ref_ids,
               const scitbx::af::shared<std::size_t> exp_ids) {
      // Expose pointers for experiment flex arrays
      const std::size_t* ptr_exp_ids = exp_ids.begin();

      std::size_t exp_it;  // ref_id;
      std::pair<const IntType*, const IntType*> refIDRange;

      unsigned int data_count = 0;
      result = 0;

      // Perform subselection of data ranges
      for (exp_it = 0; exp_it < exp_ids.size(); ++exp_it) {
        refIDRange =
          std::equal_range(ref_ids.begin(), ref_ids.end(), ptr_exp_ids[exp_it]);
        data_count = refIDRange.second - refIDRange.first;
        result += data_count;
      }
    }
    int result;
  };

  void export_pg_surpl_iter() {
    typedef return_value_policy<return_by_value> rbv;
    // templated constructor for int and size_t flex arrays of id
    class_<dials::refinement::boost_python::pg_surpl_iter>("pg_surpl_iter", no_init)
      .def(init<scitbx::af::shared<std::size_t>,
                const scitbx::af::shared<std::size_t>,
                const boost::python::list,
                const scitbx::af::shared<std::size_t>,
                const int>((boost::python::arg("ref_ids"),
                            boost::python::arg("panel_id"),
                            boost::python::arg("pnl_ids"),
                            boost::python::arg("exp_ids"),
                            boost::python::arg("cutoff"))))
      .def(init<scitbx::af::shared<int>,
                const scitbx::af::shared<std::size_t>,
                const boost::python::list,
                const scitbx::af::shared<std::size_t>,
                const int>((boost::python::arg("ref_ids"),
                            boost::python::arg("panel_id"),
                            boost::python::arg("pnl_ids"),
                            boost::python::arg("exp_ids"),
                            boost::python::arg("cutoff"))))
      .add_property(
        "result",
        make_getter(&dials::refinement::boost_python::pg_surpl_iter::result, rbv()));
  }

  void export_uc_surpl_iter() {
    typedef return_value_policy<return_by_value> rbv;
    class_<dials::refinement::boost_python::uc_surpl_iter>("uc_surpl_iter", no_init)
      .def(init<const scitbx::af::shared<std::size_t>,
                const scitbx::af::shared<cctbx::miller::index<> >,
                const scitbx::af::shared<std::size_t>,
                const scitbx::af::shared<mat3<double> > >(
        (boost::python::arg("ref_ids"),
         boost::python::arg("ref_miller"),
         boost::python::arg("exp_ids"),
         boost::python::arg("F_dbdp"))))
      .def(init<const scitbx::af::shared<int>,
                const scitbx::af::shared<cctbx::miller::index<> >,
                const scitbx::af::shared<std::size_t>,
                const scitbx::af::shared<mat3<double> > >(
        (boost::python::arg("ref_ids"),
         boost::python::arg("ref_miller"),
         boost::python::arg("exp_ids"),
         boost::python::arg("F_dbdp"))))
      .add_property(
        "result",
        make_getter(&dials::refinement::boost_python::uc_surpl_iter::result, rbv()));
  }

  void export_surpl_iter() {
    typedef return_value_policy<return_by_value> rbv;
    class_<dials::refinement::boost_python::surpl_iter>("surpl_iter", no_init)
      .def(
        init<scitbx::af::shared<std::size_t>, const scitbx::af::shared<std::size_t> >(
          (boost::python::arg("ref_ids"), boost::python::arg("exp_ids"))))
      .def(init<scitbx::af::shared<int>, const scitbx::af::shared<std::size_t> >(
        (boost::python::arg("ref_ids"), boost::python::arg("exp_ids"))))
      .add_property(
        "result",
        make_getter(&dials::refinement::boost_python::surpl_iter::result, rbv()));
  }
}}}  // namespace dials::refinement::boost_python
