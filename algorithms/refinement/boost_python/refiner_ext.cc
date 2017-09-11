#include <dials/array_family/flex_table.h>
#include <dials/array_family/reflection_table.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

#include <scitbx/matrix/transpose_multiply.h>

using namespace boost::python;
using namespace dials::af;

namespace dials { namespace refinement { namespace boost_python {

  /*
    Class used to replace functionality from panel_gp_nparam_minus_nref(p, pnl_ids, group, reflections, verbose=False) function in refiner.py
  */
  struct pgnmn_iter{
    pgnmn_iter( const reflection_table ref_table,
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

      const scitbx::af::shared<int> &ref_ids = ref_table.get<int>("id");
      //const int* ptr_ref_ids = ref_ids.begin();

      const scitbx::af::shared<size_t> &panel_id = ref_table.get<size_t>("panel");
      const std::size_t* ptr_panel_id = panel_id.begin();

      //Predefine loop variables to allow for OpenMP use;
      std::size_t exp_it, pnl; //ref_id; 

      std::pair< const int * ,        const int*         > refIDRange;
      std::pair< const std::size_t *, const std::size_t* > panelIDRange;
      unsigned int data_offset = 0;
      unsigned int data_range = 0;
      
      //Search for the beginning and end indices for the index values and track their pointers. Use these pointers to narow the next search margin and count the range of elements that match. These values indicate the count and can be added to a running total.
      for(exp_it = 0; exp_it < exp_ids.size(); ++exp_it){ //Iterate through experiment IDs
        //Get the pointer indices to the beginning and end of the values with the requested experimental ID. Performs search in ~2logN
        refIDRange = std::equal_range(ref_ids.begin(), ref_ids.end(), ptr_exp_ids[exp_it]);

        //Calculate pointer offset for start and end of the data to examine.
        data_range = refIDRange.second - refIDRange.first;
        data_offset = refIDRange.first - ref_ids.begin();
        
        //Taking the previous pointer positions, and perfoming a search within this range for the matching panel IDs, we then counti the pointer range as all valid entries and add to the accumlator.
        for(pnl = 0; pnl < pnl_ids_flex.size(); ++pnl){ //Iterate through panel IDs
          panelIDRange = std::equal_range( ptr_panel_id + data_offset, 
                                           ptr_panel_id + data_offset + data_range, 
                                           pnl_ids_flex[pnl] );
          count += (unsigned int)(panelIDRange.second - panelIDRange.first);
        }
      }
      result = count - cutoff;
    }
    int result;
  };

  /*
    Class used to replace functionality from unit_cell_nparam_minus_nref(p, reflections) function in refiner.py
  */
  struct ucnmn_iter{
    ucnmn_iter( const reflection_table ref_table,
                const scitbx::af::shared<std::size_t> exp_ids, 
                const scitbx::af::shared<mat3<double> > F_dbdp
              )
    {
      //Expose pointers for experiment flex arrays
      const std::size_t* ptr_exp_ids = exp_ids.begin();
      const scitbx::af::shared<int> &ref_ids = ref_table.get<int>("id");

      const scitbx::af::shared<cctbx::miller::index<> > &ref_miller = ref_table.get<cctbx::miller::index<> >("miller_index");
      
      std::size_t exp_it; //ref_id; 
      
      std::pair< const int * ,const int * > refIDRange;

      unsigned int data_range = 0;
      unsigned int data_offset = 0;
      
      //Perform subselection of data ranges
      for( exp_it = 0; exp_it < exp_ids.size(); ++exp_it ){
        refIDRange = std::equal_range( ref_ids.begin(),
                                       ref_ids.end(),
                                       ptr_exp_ids[exp_it] );
        data_range = refIDRange.second - refIDRange.first;
        data_offset = refIDRange.first - ref_ids.begin();
      }

      //Pre-initialise loop variables and array for storing norms
      int ii, jj;
      scitbx::af::shared< double > mul_norm;

      //Calculate L-2 norm of Matrix*Vec operation
      for( ii=0; ii< data_range; ++ii){
        vec3<double> v3(ref_miller[ii+data_offset]);
        for( jj=0; jj < F_dbdp.size(); ++jj){
          mul_norm.push_back( (F_dbdp[jj]*v3).length() );
        }
      }

      //Count the number of elements greater than 0 from the L-2 norms.
      scitbx::af::shared<int> nref_each_param;
      for (ii=0; ii < F_dbdp.size(); ++ii){
        nref_each_param.push_back(
          std::count_if (mul_norm.begin() + ii*data_range, 
                         mul_norm.begin() + (ii+1)*data_range, 
                         gtZero) 
          );
      }

      //Return the value of the smallest element
      scitbx::af::shared<int>::iterator it = std::min_element(nref_each_param.begin(),nref_each_param.end());
      result = *(nref_each_param.begin()+std::distance(std::begin(nref_each_param), it));
    }
    int result;
    static bool gtZero (double x) { return ( (x > 0.0) ==1 ); } //Used for comparison operation
  };

  void export_pgnmn_iter()
  {
    typedef return_value_policy<return_by_value> rbv;

    class_<dials::refinement::boost_python::pgnmn_iter>
    ("pgnmn_iter", init<reflection_table,
      const boost::python::list,
      const scitbx::af::shared<std::size_t>,
      const int>
      ((
        boost::python::arg("ref_table"),
        boost::python::arg("pnl_ids"),
        boost::python::arg("exp_ids"),
        boost::python::arg("cutoff")
      ))
    ).add_property("result",make_getter(&dials::refinement::boost_python::pgnmn_iter::result, rbv()))
  ;
  }

  void export_ucnmn_iter()
  {
    typedef return_value_policy<return_by_value> rbv;

    class_<dials::refinement::boost_python::ucnmn_iter>
    ("ucnmn_iter", init<
      const reflection_table,
      const scitbx::af::shared<std::size_t>, 
      const scitbx::af::shared<mat3<double> > >
      ((
        boost::python::arg("ref_table"),
        boost::python::arg("exp_ids"),
        boost::python::arg("dbdp"),
        boost::python::arg("len_dbdp")
      ))
    ).add_property("result",make_getter(&dials::refinement::boost_python::ucnmn_iter::result, rbv()))
    ;
  }
}}} // namespace dials::refinement::boost_python
