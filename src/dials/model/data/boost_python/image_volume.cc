/*
 * image_volume.cc
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
#include <scitbx/array_family/flex_types.h>
#include <dials/array_family/reflection_table.h>
#include <dials/model/data/image_volume.h>
#include <dials/model/data/mask_code.h>
#include <dials/config.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using af::BackgroundIncludesBadPixels;
  using af::ForegroundIncludesBadPixels;

  template <typename FloatType>
  void MultiPanelImageVolume_update_reflection_info(
    MultiPanelImageVolume<FloatType> image_volume,
    af::reflection_table reflections) {
    // Check the input
    DIALS_ASSERT(reflections.contains("bbox"));
    DIALS_ASSERT(reflections.contains("panel"));

    // Get some reflection table columns
    af::const_ref<int6> bbox = reflections["bbox"];
    af::const_ref<std::size_t> panel = reflections["panel"];

    // Get the flags array
    af::ref<std::size_t> flags = reflections["flags"];

    // Set some information about number of pixels
    af::ref<std::size_t> num_valid = reflections["num_pixels.valid"];
    af::ref<std::size_t> num_bg = reflections["num_pixels.background"];
    // af::ref<std::size_t> num_bg_used = reflections["num_pixels.background_used"];
    af::ref<std::size_t> num_fg = reflections["num_pixels.foreground"];
    af::ref<std::size_t> num_ol = reflections["num_pixels.overlapped"];

    // Set some background information
    af::ref<double> bg_mean = reflections["background.mean"];
    af::ref<double> bg_disp = reflections["background.dispersion"];

    // The mask codes
    uint8_t mask_code1 = Valid;
    uint8_t mask_code2 = Valid | Background;
    uint8_t mask_code3 = Valid | Background | BackgroundUsed;
    uint8_t mask_code4 = Valid | Foreground;
    uint8_t mask_code5 = Overlapped;

    // Loop through each reflection
    for (std::size_t i = 0; i < panel.size(); ++i) {
      // Get the image volume
      ImageVolume<FloatType> v = image_volume.get(panel[i]);
      DIALS_ASSERT(v.is_consistent());

      // Trim the bounding box
      int6 b = v.trim_bbox(bbox[i]);

      // Get the data arrays
      af::versa<FloatType, af::c_grid<3>> data = v.extract_data(b);
      af::versa<FloatType, af::c_grid<3>> bgrd = v.extract_background(b);
      af::versa<uint8_t, af::c_grid<3>> mask = v.extract_mask(b, i);

      // Compute numbers of pixels
      std::size_t num1 = 0;
      std::size_t num2 = 0;
      std::size_t num3 = 0;
      std::size_t num4 = 0;
      std::size_t num5 = 0;
      std::size_t num6 = 0;
      std::size_t num7 = 0;
      for (std::size_t j = 0; j < mask.size(); ++j) {
        if ((mask[j] & mask_code1) == mask_code1) num1++;
        if ((mask[j] & mask_code2) == mask_code2) num2++;
        if ((mask[j] & mask_code3) == mask_code3) num3++;
        if ((mask[j] & mask_code4) == mask_code4) num4++;
        if ((mask[j] & mask_code5) == mask_code5) num5++;
        if ((mask[j] & Background) && !(mask[j] & Valid)) num6++;
        if ((mask[j] & Foreground) && !(mask[j] & Valid)) num7++;
      }
      num_valid[i] = num1;
      num_bg[i] = num2;
      // num_bg_used[i] = num3;
      num_fg[i] = num4;
      num_ol[i] = num5;

      // Set some flags
      if (num6 > 0) {
        flags[i] |= BackgroundIncludesBadPixels;
      } else {
        flags[i] &= ~BackgroundIncludesBadPixels;
      }
      if (num7 > 0) {
        flags[i] |= ForegroundIncludesBadPixels;
      } else {
        flags[i] &= ~ForegroundIncludesBadPixels;
      }

      // Compute the mean (modelled) background and dispersion in background
      double fg_sum = 0.0;
      double fg_cnt = 0.0;
      double bg_sum = 0.0;
      double bg_cnt = 0.0;
      double bg_sum_sq = 0.0;
      for (std::size_t j = 0; j < mask.size(); ++j) {
        if (mask[j] & Foreground) {
          fg_sum += bgrd[j];
          fg_cnt += 1;
        }
        if ((mask[j] | mask_code2) == mask_code2) {
          bg_sum += data[j];
          bg_sum_sq += data[j] * data[j];
          bg_cnt += 1;
        }
      }
      if (fg_cnt > 0) {
        bg_mean[i] = fg_sum / fg_cnt;
      } else {
        bg_mean[i] = 0.0;
      }
      if (bg_cnt > 0 && bg_sum > 0) {
        bg_disp[i] = bg_sum_sq / bg_sum - bg_sum / bg_cnt;
      } else {
        bg_disp[i] = 0.0;
      }
    }
  }

  template <typename FloatType>
  void image_volume_wrapper(const char* name) {
    typedef ImageVolume<FloatType> Class;

    struct IVPickleSuite : pickle_suite {
      // The mask is pickled as a flex_grid array, as the c_grid<3> conversions for
      // uint8_t are only registered once dials.array_family.flex has been imported
      typedef af::versa<uint8_t, af::flex_grid<>> flex_mask_type;

      static tuple getinitargs(const Class& obj) {
        af::c_grid<3> grid = obj.accessor();
        return make_tuple(obj.frame0(), obj.frame1(), grid[1], grid[2]);
      }

      static tuple getstate(const Class& obj) {
        af::versa<uint8_t, af::c_grid<3>> mask = obj.mask();
        af::versa<Label, af::c_grid<3>> label = obj.label1();
        af::versa<int, af::c_grid<3>> label_first(label.accessor());
        af::versa<int, af::c_grid<3>> label_second(label.accessor());
        for (std::size_t i = 0; i < label.size(); ++i) {
          label_first[i] = label[i].first;
          label_second[i] = label[i].second;
        }
        return make_tuple(obj.data(),
                          obj.background(),
                          flex_mask_type(mask, mask.accessor().as_flex_grid()),
                          label_first,
                          label_second);
      }

      static void setstate(Class& obj, tuple state) {
        DIALS_ASSERT(len(state) == 5);
        const flex_mask_type& mask = extract<const flex_mask_type&>(state[2]);
        obj.set_state(extract<af::const_ref<FloatType, af::c_grid<3>>>(state[0]),
                      extract<af::const_ref<FloatType, af::c_grid<3>>>(state[1]),
                      af::const_ref<uint8_t, af::c_grid<3>>(
                        mask.begin(), af::c_grid<3>(mask.accessor())),
                      extract<af::const_ref<int, af::c_grid<3>>>(state[3]),
                      extract<af::const_ref<int, af::c_grid<3>>>(state[4]));
      }
    };

    class_<Class>(name, no_init)
      .def(init<int, int, std::size_t, std::size_t>())
      .def("frame0", &Class::frame0)
      .def("frame1", &Class::frame1)
      .def("accessor", &Class::accessor)
      .def("data", &Class::data)
      .def("background", &Class::background)
      .def("mask", &Class::mask)
      .def("set_image", &Class::template set_image<int>)
      .def("set_image", &Class::template set_image<double>)
      .def("is_consistent", &Class::is_consistent)
      .def("extract_data", &Class::extract_data)
      .def("extract_background", &Class::extract_background)
      .def("extract_mask", &Class::extract_mask)
      .def_pickle(IVPickleSuite());
  }

  template <typename FloatType>
  void multi_panel_image_volume_wrapper(const char* name) {
    typedef MultiPanelImageVolume<FloatType> Class;

    struct MPIVPickleSuite : pickle_suite {
      static tuple getstate(const Class& obj) {
        list volumes;
        for (std::size_t i = 0; i < obj.size(); ++i) {
          volumes.append(obj.get(i));
        }
        return tuple(volumes);
      }

      static void setstate(Class& obj, tuple state) {
        DIALS_ASSERT(obj.size() == 0);
        for (std::size_t i = 0; i < len(state); ++i) {
          obj.add(extract<const ImageVolume<FloatType>&>(state[i]));
        }
      }
    };

    class_<Class>(name)
      .def("frame0", &Class::frame0)
      .def("frame1", &Class::frame1)
      .def("add", &Class::add)
      .def("get", &Class::get)
      .def("set_image", &Class::template set_image<int>)
      .def("set_image", &Class::template set_image<double>)
      .def("update_reflection_info",
           &MultiPanelImageVolume_update_reflection_info<FloatType>)
      .def("__len__", &Class::size)
      .def_pickle(MPIVPickleSuite());
  }

  void export_image_volume() {
    typedef ImageVolume<>::float_type FloatType;

    image_volume_wrapper<FloatType>("ImageVolume");
    multi_panel_image_volume_wrapper<FloatType>("MultiPanelImageVolume");
  }

}}}  // namespace dials::model::boost_python
