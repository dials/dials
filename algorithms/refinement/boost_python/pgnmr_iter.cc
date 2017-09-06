#include <dials/array_family/flex_table.h>
#include <dials/array_family/reflection_table.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;
using namespace dials::af;

namespace dials { namespace refinement { namespace boost_python {

  struct pgnmr_iter{
    pgnmr_iter( reflection_table ref_table,
                const boost::python::list pnl_ids,
                const scitbx::af::shared<std::size_t> exp_ids,
                const int cutoff){
      //Perform python list conversion once; improves read access
      scitbx::af::shared<std::size_t> pnl_ids_flex;
      for (int ii=0; ii < len(pnl_ids); ++ii)
        pnl_ids_flex.push_back(boost::python::extract<std::size_t>(pnl_ids[ii]));

      std::size_t count = 0;

      //Work with raw pointers instead of overloadable objects
      const std::size_t* ptr_exp_ids = exp_ids.begin();

      scitbx::af::shared<int> ref_ids = ref_table["id"];
      const int* ptr_ref_ids = ref_ids.begin();

      scitbx::af::shared<size_t> panel_id = ref_table["panel"];
      const std::size_t* ptr_panel_id = panel_id.begin();

      //Predefine loop variables to allow for OpenMP use;
      std::size_t exp_it, ref_id, pnl;

      //Will give performance benefit if OpenMP available;
      #pragma omp parallel for private(ref_id,pnl) schedule(dynamic) reduction(+:count)
      for ( exp_it = 0; exp_it < exp_ids.size(); ++exp_it ){ //Iterate through experiment IDs
        for( ref_id = 0; ref_id < ref_ids.size(); ++ref_id ){ //Iterate through reflection table IDs
          if ( ptr_ref_ids[ref_id] == ptr_exp_ids[exp_it] ){ //Compare reflection table ID with given experiment IDs;
            for (pnl = 0; pnl < pnl_ids_flex.size(); ++pnl){ //Iterate through panel IDs
              if( pnl_ids_flex[pnl] == ptr_panel_id[ref_id]){ //Compare panel IDs for the selected experiments against the given list
                ++count; //Increment if the above matches are true.
              }
            }
          }
        }
      }
      result = count - cutoff;
    }
    int result;
  };
  void export_pgnmr_iter()
  {
    typedef return_value_policy<return_by_value> rbv;

    class_<dials::refinement::boost_python::pgnmr_iter>
    ("pgnmr_iter", init<reflection_table,
      const boost::python::list,
      const scitbx::af::shared<std::size_t>,
      const int>
      ((
        boost::python::arg("ref_table"),
        boost::python::arg("pnl_ids"),
        boost::python::arg("exp_ids"),
        boost::python::arg("cutoff")
      ))
    ).add_property("result",make_getter(&dials::refinement::boost_python::pgnmr_iter::result, rbv()))
  ;
  }
}}} // namespace dials::refinement::boost_python
