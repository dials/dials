/*
 * flex_shoebox.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/unordered_map.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <cmath>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/partial_shoebox.h>
#include <dials/model/data/shoebox.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using af::int2;
  using af::int6;
  using af::small;
  using dials::model::PartialShoebox;
  using dials::model::Shoebox;

  /**
   * Check if the arrays are consistent
   */
  static
  shared<bool> is_consistent(const const_ref<PartialShoebox> &a) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_consistent();
    }
    return result;
  }

  /**
   * Get the bounding boxes
   */
  static
  shared<int6> bounding_boxes(const const_ref<PartialShoebox> &a) {
    shared<int6> result(a.size(), af::init_functor_null<int6>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].bbox;
    }
    return result;
  }
  
  /**
   * Get the panel numbers
   */
  static
  shared<std::size_t> panels(const const_ref<PartialShoebox> &a) {
    shared<std::size_t> result(a.size(), af::init_functor_null<std::size_t>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].panel;
    }
    return result;
  }
  
  /**
   * Get the zranges
   */
  static
  shared<int2> zranges(const const_ref<PartialShoebox> &a) {
    shared<int2> result(a.size(), af::init_functor_null<int2>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].zrange;
    }
    return result;
  }
  
  /**
   * Struct for sorting shoeboxes
   */
  struct CompareMinZ {
    const_ref<PartialShoebox> x_;
    CompareMinZ(const const_ref<PartialShoebox> &x) : x_(x) {}
    bool operator()(int a, int b) const {
      return x_[a].zrange[0] < x_[b].zrange[0];
    }
  };
   
  /**
   * Try to merge the partial shoeboxes
   */
  template <typename FloatType>
  Shoebox<FloatType> merge(const const_ref<PartialShoebox> &a, int2 srange) {
  
    // Sort the array by minimum z
    DIALS_ASSERT(a.size() > 0);
    DIALS_ASSERT(srange[1] > srange[0]);
    shared<int> index(a.size());
    for (std::size_t i = 0; i < index.size(); ++i) index[i] = i;
    std::sort(index.begin(), index.end(), CompareMinZ(a));

    // Get the bbox and panel
    int6 bbox = a[0].bbox;
    std::size_t panel = a[0].panel;
        
    // Loop through and check ranges are valid and within the bbox range
    DIALS_ASSERT(a[index[0]].is_consistent());
    DIALS_ASSERT(a[index[0]].zrange[0] == std::max(bbox[4], srange[0]));
    for (std::size_t i = 1; i < index.size(); ++i) {
      DIALS_ASSERT(a[index[i]].is_consistent());
      DIALS_ASSERT(a[index[i]].bbox.all_eq(bbox));
      DIALS_ASSERT(a[index[i]].panel == panel);
      int2 z0 = a[index[i-1]].zrange;
      int2 z1 = a[index[i]].zrange;
      DIALS_ASSERT(z1[0] == z0[1]);
    }
    DIALS_ASSERT(a[index.back()].zrange[1] <= std::min(bbox[5], srange[1]));

    // Copy all the data
    Shoebox<FloatType> shoebox(panel, bbox);
    shoebox.allocate();
    std::size_t xysize = shoebox.xsize() * shoebox.ysize();
    for (std::size_t i = 0; i < a.size(); ++i) {
      std::size_t offset = (a[i].zrange[0] - bbox[4]) * xysize;
      for (std::size_t j = 0; j < a[i].data.size(); ++j) {
        shoebox.data[offset + j] = (FloatType)a[i].data[j];
      }
    }
    
    // Return the new shoebox
    return shoebox;
  }
  
  /**
   * Merge all shoeboxes and return indices and shoeboxes
   */
  template <typename FloatType>
  std::pair< af::shared<std::size_t>, af::shared< Shoebox<FloatType> > >
  merge_all(const const_ref<PartialShoebox> &shoeboxes, 
      const const_ref<std::size_t> &indices, int2 srange) {

    // Some useful typedefs
    typedef boost::unordered_map<std::size_t, af::shared<PartialShoebox> > map_type;
    typedef typename map_type::iterator iterator;
    typedef Shoebox<FloatType> shoebox_type;

    // Construct a map of the shoeboxes to merge
    DIALS_ASSERT(indices.size() == shoeboxes.size());
    map_type tomerge(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      tomerge[indices[i]].push_back(shoeboxes[i]);
    }
    
    // Merge all the shoeboxes
    af::shared<std::size_t> cindices;
    af::shared<shoebox_type> cshoeboxes;
    for (iterator it = tomerge.begin(); it != tomerge.end(); ++it) {
      cindices.push_back(it->first);
      cshoeboxes.push_back(merge<FloatType>(it->second.const_ref(), srange));
    }
    return std::make_pair(cindices, cshoeboxes);
  }
   
  /**
   * A class to convert the shoebox class to a string for pickling
   */
  struct partial_shoebox_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    /** Initialise with the version for checking */
    partial_shoebox_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    partial_shoebox_to_string& operator<<(const PartialShoebox &val) {
      *this << val.panel
            << val.bbox[0]
            << val.bbox[1]
            << val.bbox[2]
            << val.bbox[3]
            << val.bbox[4]
            << val.bbox[5]
            << val.zrange[0]
            << val.zrange[1];
            
      profile_to_string(val.data);

      return *this;
    }
    
    /** Convert a profile to string */
    template <typename ProfileType>
    void profile_to_string(const ProfileType &p) {
      *this << p.accessor().size();
      for (std::size_t i = 0; i < p.accessor().size(); ++i) {
        *this << p.accessor()[i];
      }
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this << p[i];
      }
    }
  };

  /** 
   * A class to convert a string to a shoebox for unpickling
   */
  struct partial_shoebox_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    /** Initialise the class with the string. Get the version and check */
    partial_shoebox_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    partial_shoebox_from_string& operator>>(PartialShoebox &val) {
      *this >> val.panel
            >> val.bbox[0]
            >> val.bbox[1]
            >> val.bbox[2]
            >> val.bbox[3]
            >> val.bbox[4]
            >> val.bbox[5]
            >> val.zrange[0]
            >> val.zrange[1];

      val.data = profile_from_string< versa<int, c_grid<3> > >();

      return *this;
    }

    /** Get a profile from a string */
    template <typename ProfileType>
    ProfileType profile_from_string() {
      typename ProfileType::accessor_type accessor;
      typename ProfileType::size_type n_dim;
      *this >> n_dim;
      DIALS_ASSERT(n_dim == accessor.size());
      for (std::size_t i = 0; i < n_dim; ++i) {
        *this >> accessor[i];
      }
      ProfileType p = ProfileType(accessor);
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this >> p[i];
      }
      return p;
    }

    unsigned int version;
  };

  scitbx::af::boost_python::flex_wrapper<
    PartialShoebox,
    return_internal_reference<> >::class_f_t
  flex_partial_shoebox_wrapper(const char *name)
  {
    return scitbx::af::boost_python::flex_wrapper <
      PartialShoebox, return_internal_reference<> >::plain(name)
        .def("is_consistent", &is_consistent)
        .def("panels", &panels)
        .def("bounding_boxes", &bounding_boxes)
        .def("zranges", &zranges)
        .def("merge", &merge<ProfileFloatType>)
        .def("merge_all", &merge_all<ProfileFloatType>)
        .def_pickle(flex_pickle_double_buffered<PartialShoebox, 
          partial_shoebox_to_string, 
          partial_shoebox_from_string>());
  }

  void export_flex_partial_shoebox() {
    flex_partial_shoebox_wrapper("partial_shoebox");

    boost_adaptbx::std_pair_conversions::to_and_from_tuple<
      af::shared<std::size_t>, 
      af::shared<Shoebox<ProfileFloatType> > >();

  }

}}} // namespace dials::af::boost_python
