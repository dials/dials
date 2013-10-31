/*
 * flex_pixel_list.cc
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
#include <vector>
#include <algorithm>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/pixel_list.h>
#include <dials/error.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using scitbx::vec2;
  using scitbx::vec3;
  using dials::model::PixelList;
  
  struct PixelListLessPred {
    PixelListLessPred(const af::const_ref<PixelList> a)
      : pl_(a) {}
     
    template <typename T> 
    bool operator()(T a, T b) const {
      return pl_[a].first_frame() < pl_[b].last_frame();
    }
    
    af::const_ref<PixelList> pl_;
  };
  
  PixelList pixel_list_merge(const af::const_ref<PixelList> &a) {
    
    // Sort by frames  
    std::vector<int> index(a.size());
    for (std::size_t i = 0; i < index.size(); ++i) {
      index[i] = i;
    }
    std::sort(index.begin(), index.end(), PixelListLessPred(a));

    std::size_t num = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      num += a[i].num_pixels();
    }
    
    af::shared< vec3<int> > coords((af::reserve(num)));
    af::shared<int> values((af::reserve(num)));
    
    // Check the frames are contiguous
    int2 fr1(0, a[index[0]].first_frame());
    int2 size = a[index[0]].size();
    for (std::size_t i = 1; i < index.size(); ++i) {
      const PixelList &pl = a[index[i]];
      int2 fr2 = pl.frame_range();
      DIALS_ASSERT(fr2[0] == fr1[1]);
      DIALS_ASSERT(pl.size().all_eq(size));
      fr1 = fr2;
      af::shared<int> v = pl.values();
      af::shared< vec3<int> > c = pl.coords();        
      std::copy(v.begin(), v.end(), std::back_inserter(values));
      std::copy(c.begin(), c.end(), std::back_inserter(coords));
    }
    
    int2 frame_range(a[index.front()].first_frame(), a[index.back()].last_frame());
    
    return PixelList(size, frame_range, values, coords); 
  }
  
  void export_flex_pixel_list()
  {
    scitbx::af::boost_python::flex_wrapper <PixelList, 
        return_internal_reference<> >::plain("pixel_list")
      .def("merge", &pixel_list_merge);
  }

}}} // namespace dials::af::boost_python
