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
#include <cmath>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/pixel_list.h>
#include <dials/model/data/observation.h>
#include <dials/algorithms/image/connected_components/connected_components.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using af::int2;
  using af::int6;
  using af::small;
  using scitbx::vec3;
  using dials::model::Shoebox;
  using dials::model::Centroid;
  using dials::model::Intensity;
  using dials::model::Observation;
  using dials::model::PixelList;
  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::model::BackgroundUsed;
  using dials::algorithms::LabelImageStack;
  using dials::algorithms::LabelPixels;

  /**
   * Construct from an array of panels and bounding boxes.
   */
  template <typename FloatType>
  typename af::flex< Shoebox<FloatType> >::type* from_panel_and_bbox(
      const af::const_ref<std::size_t> panel,
      const af::const_ref<int6> bbox,
      bool allocate,
      bool flatten) {
    DIALS_ASSERT(panel.size() == bbox.size());
    af::shared< Shoebox<FloatType> > result(panel.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = Shoebox<FloatType>(panel[i], bbox[i], flatten);
      if (allocate) {
        result[i].allocate();
      }
    }
    return new typename af::flex< Shoebox<FloatType> >::type(
      result, af::flex_grid<>(result.size()));
  }

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <typename FloatType>
  typename af::flex< Shoebox<FloatType> >::type* from_pixel_list(
      const PixelList &pixel, std::size_t panel,
      std::size_t zstart, bool twod) {

    // Get the stuff from the label struct
    af::shared<int> labels = twod ? pixel.labels_2d() : pixel.labels_3d();
    af::shared<double> values = pixel.values();
    af::shared< vec3<int> > coords = pixel.coords();

    // Get the number of labels and allocate the array
    std::size_t num = af::max(labels.const_ref()) + 1;
    af::shared< Shoebox<FloatType> > result(num, Shoebox<FloatType>());

    // Initialise the bboxes
    int xsize = pixel.size()[1];
    int ysize = pixel.size()[0];
    int2 minmaxz = pixel.frame_range();
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].bbox[0] = xsize; result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize; result[i].bbox[3] = 0;
      result[i].bbox[4] = minmaxz[1]; result[i].bbox[5] = minmaxz[0];
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] <  result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] <  result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] <  result[l].bbox[4]) result[l].bbox[4] = c[0];
      if (c[0] >= result[l].bbox[5]) result[l].bbox[5] = c[0] + 1;
    }

    // Allocate all the arrays
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].allocate();
    }

    // Set all the mask and data points
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      FloatType v = values[i];
      vec3<int> c = coords[i];
      int ii = c[2] - result[l].bbox[0];
      int jj = c[1] - result[l].bbox[2];
      int kk = c[0] - result[l].bbox[4];
      DIALS_ASSERT(ii >= 0 && jj >= 0 && kk >= 0);
      DIALS_ASSERT(ii < result[l].xsize());
      DIALS_ASSERT(jj < result[l].ysize());
      DIALS_ASSERT(kk < result[l].zsize());
      result[l].data(kk,jj,ii) = v;
      result[l].mask(kk,jj,ii) = Valid | Foreground;
    }

    // Shift bbox z start position
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].bbox[4] += zstart;
      result[i].bbox[5] += zstart;
    }

    // Return the array
    return new typename af::flex< Shoebox<FloatType> >::type(
      result, af::flex_grid<>(num));
  }

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <std::size_t DIM, typename FloatType>
  typename af::flex< Shoebox<FloatType> >::type* from_labels(
      const LabelImageStack<DIM> &label,
      std::size_t panel, std::size_t zstart) {

    // Get the stuff from the label struct
    af::shared<int> labels = label.labels();
    af::shared<int> values = label.values();
    af::shared< vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = af::max(labels.const_ref()) + 1;
    af::shared< Shoebox<FloatType> > result(num, Shoebox<FloatType>());

    // Initialise the bboxes
    int xsize = label.size()[1];
    int ysize = label.size()[0];
    int zsize = label.num_images();
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].bbox[0] = xsize; result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize; result[i].bbox[3] = 0;
      result[i].bbox[4] = zsize; result[i].bbox[5] = 0;
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] <  result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] <  result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] <  result[l].bbox[4]) result[l].bbox[4] = c[0];
      if (c[0] >= result[l].bbox[5]) result[l].bbox[5] = c[0] + 1;
    }

    // Allocate all the arrays
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].allocate();
    }

    // Set all the mask and data points
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      FloatType v = values[i];
      vec3<int> c = coords[i];
      int ii = c[2] - result[l].bbox[0];
      int jj = c[1] - result[l].bbox[2];
      int kk = c[0] - result[l].bbox[4];
      DIALS_ASSERT(ii >= 0 && jj >= 0 && kk >= 0);
      DIALS_ASSERT(ii < result[l].xsize());
      DIALS_ASSERT(jj < result[l].ysize());
      DIALS_ASSERT(kk < result[l].zsize());
      result[l].data(kk,jj,ii) = v;
      result[l].mask(kk,jj,ii) = Valid | Foreground;
    }

    // Shift bbox z start position
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].bbox[4] += zstart;
      result[i].bbox[5] += zstart;
    }

    // Return the array
    return new typename af::flex< Shoebox<FloatType> >::type(
      result, af::flex_grid<>(num));
  }

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <typename FloatType>
  typename af::flex< Shoebox<FloatType> >::type* from_pixel_labeller(
      const LabelPixels &label, std::size_t panel) {

    // Get the stuff from the label struct
    af::shared<int> labels = label.labels();
    af::shared<int> values = label.values();
    af::shared< vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = af::max(labels.const_ref()) + 1;
    af::shared< Shoebox<FloatType> > result(num, Shoebox<FloatType>());

    // Initialise the bboxes
    int xsize = label.size()[2];
    int ysize = label.size()[1];
    int zsize = label.size()[0];
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].bbox[0] = xsize; result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize; result[i].bbox[3] = 0;
      result[i].bbox[4] = zsize; result[i].bbox[5] = 0;
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] <  result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] <  result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] <  result[l].bbox[4]) result[l].bbox[4] = c[0];
      if (c[0] >= result[l].bbox[5]) result[l].bbox[5] = c[0] + 1;
    }

    // Allocate all the arrays
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].allocate();
    }

    // Set all the mask and data points
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      double v = values[i];
      vec3<int> c = coords[i];
      int ii = c[2] - result[l].bbox[0];
      int jj = c[1] - result[l].bbox[2];
      int kk = c[0] - result[l].bbox[4];
      DIALS_ASSERT(ii >= 0 && jj >= 0 && kk >= 0);
      DIALS_ASSERT(ii < result[l].xsize());
      DIALS_ASSERT(jj < result[l].ysize());
      DIALS_ASSERT(kk < result[l].zsize());
      result[l].data(kk,jj,ii) = (double)v;
      result[l].mask(kk,jj,ii) = Valid | Foreground;
    }

    // Return the array
    return new typename af::flex< Shoebox<FloatType> >::type(
      result, af::flex_grid<>(num));
  }

  /**
   * Allocate the shoeboxes
   */
  template <typename FloatType>
  void allocate(af::ref<Shoebox<FloatType> > a) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      a[i].allocate();
    }
  }

  /**
   * Allocate the shoeboxes
   */
  template <typename FloatType>
  void allocate_with_value(af::ref<Shoebox<FloatType> > a, int mask_code) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      a[i].allocate_with_value(mask_code);
    }
  }

  /**
   * Deallocate the shoeboxes
   */
  template <typename FloatType>
  void deallocate(af::ref<Shoebox<FloatType> > a) {
    for (std::size_t i = 0; i < a.size(); ++i) {
      a[i].deallocate();
    }
  }

  /**
   * Check if the arrays are consistent
   */
  template <typename FloatType>
  shared<bool> is_consistent(const const_ref< Shoebox<FloatType> > &a) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_consistent();
    }
    return result;
  }

  /**
   * Check if the bounding box has points outside the image range.
   */
  template <typename FloatType>
  shared<bool> is_bbox_within_image_volume(
      const const_ref< Shoebox<FloatType> > &a,
      int2 image_size, int2 scan_range) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_bbox_within_image_volume(image_size, scan_range);
    }
    return result;
  }

  /**
   * Check if the bounding box has points that cover bad pixels
   */
  template <typename FloatType>
  shared<bool> does_bbox_contain_bad_pixels(
      const const_ref< Shoebox<FloatType> > &a,
      const const_ref<bool, c_grid<2> > &mask) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].does_bbox_contain_bad_pixels(mask);
    }
    return result;
  }

  /**
   * Count the number of mask pixels with the given code
   */
  template <typename FloatType>
  shared<int> count_mask_values(
      const const_ref< Shoebox<FloatType> > &a, int code) {
    shared<int> result(a.size(), af::init_functor_null<int>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].count_mask_values(code);
    }
    return result;
  }

  /**
   * Get the maximum index of each shoebox
   */
  template <typename FloatType>
  shared< vec3<double> > peak_coordinates(ref< Shoebox<FloatType> > a) {
    shared< vec3<double> > result(a.size(),
      af::init_functor_null< vec3<double> >());
    for (std::size_t i = 0; i < a.size(); ++i) {
      std::size_t index = af::max_index(a[i].data.const_ref());
      af::c_grid<3> accessor = a[i].data.accessor();
      tiny<std::size_t, 3> coord = accessor.index_nd(index);
      result[i][0] = a[i].bbox[0] + coord[2] + 0.5;
      result[i][1] = a[i].bbox[2] + coord[1] + 0.5;
      result[i][2] = a[i].bbox[4] + coord[0] + 0.5;
    }
    return result;
  }

  /**
   * Get the bounding boxes
   */
  template <typename FloatType>
  shared<int6> bounding_boxes(const const_ref< Shoebox<FloatType> > &a) {
    shared<int6> result(a.size(), af::init_functor_null<int6>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].bbox;
    }
    return result;
  }

  /**
   * Get the panel numbers
   */
  template <typename FloatType>
  shared<std::size_t> panels(const const_ref< Shoebox<FloatType> > &a) {
    shared<std::size_t> result(a.size(), af::init_functor_null<std::size_t>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].panel;
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_all(const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_all();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_masked(
      const const_ref< Shoebox<FloatType> > &a, int code) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_masked(code);
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_valid(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_valid();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_foreground(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_foreground();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_strong(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_strong();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_all_minus_background(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_all_minus_background();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_masked_minus_background(
      const const_ref< Shoebox<FloatType> > &a, int code) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_masked_minus_background(code);
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_valid_minus_background(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_valid_minus_background();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_foreground_minus_background(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_foreground_minus_background();
    }
    return result;
  }

  /**
   * Get a list of centroid
   */
  template <typename FloatType>
  af::shared<Centroid> centroid_strong_minus_background(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Centroid> result(a.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].centroid_strong_minus_background();
    }
    return result;
  }

  /**
   * Get a list of intensities
   */
  template <typename FloatType>
  af::shared<Intensity> summed_intensity(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<Intensity> result(a.size(), Intensity());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].summed_intensity();
    }
    return result;
  }

  /**
   * Get the mean background.
   */
  template <typename FloatType>
  af::shared<double> mean_background(
      const const_ref< Shoebox<FloatType> > &a) {
    af::shared<double> result(a.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      af::versa<FloatType, af::c_grid<3> > data = a[i].data;
      af::versa<int, af::c_grid<3> > mask = a[i].mask;
      double mean = 0.0;
      std::size_t count = 0;
      for (std::size_t j = 0; j < data.size(); ++j) {
        if (mask[j] & BackgroundUsed) {
          mean += data[j];
          count += 1;
        }
      }
      if (count > 0) {
        mean /= count;
      } else {
        mean = 0;
      }
      result[i] = mean;
    }
    return result;
  }

  /**
   * Flatten the shoeboxes
   */
  template <typename FloatType>
  void flatten(ref< Shoebox<FloatType> > self) {
    for (std::size_t i = 0; i < self.size(); ++i) {
      self[i].flatten();
    }
  }

  /**
   * A class to convert the shoebox class to a string for pickling
   */
  template <typename FloatType>
  struct shoebox_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    typedef Shoebox<FloatType> shoebox_type;

    /** Initialise with the version for checking */
    shoebox_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    shoebox_to_string& operator<<(const shoebox_type &val) {
      *this << val.panel
            << val.bbox[0]
            << val.bbox[1]
            << val.bbox[2]
            << val.bbox[3]
            << val.bbox[4]
            << val.bbox[5];

      profile_to_string(val.data);
      profile_to_string(val.mask);
      profile_to_string(val.background);

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
  template <typename FloatType>
  struct shoebox_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    typedef Shoebox<FloatType> shoebox_type;

    /** Initialise the class with the string. Get the version and check */
    shoebox_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    shoebox_from_string& operator>>(shoebox_type &val) {
      *this >> val.panel
            >> val.bbox[0]
            >> val.bbox[1]
            >> val.bbox[2]
            >> val.bbox[3]
            >> val.bbox[4]
            >> val.bbox[5];

      val.data = profile_from_string< versa<FloatType, c_grid<3> > >();
      val.mask = profile_from_string< versa<int, c_grid<3> > >();
      val.background = profile_from_string< versa<FloatType, c_grid<3> > >();

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

  template <typename FloatType>
  typename scitbx::af::boost_python::flex_wrapper<
    Shoebox<FloatType>,
    return_internal_reference<> >::class_f_t
  flex_shoebox_wrapper(const char *name)
  {
    typedef Shoebox<FloatType> shoebox_type;

    return scitbx::af::boost_python::flex_wrapper <
      shoebox_type, return_internal_reference<> >::plain(name)
        .def("__init__", make_constructor(
          from_pixel_list<FloatType>,
          default_call_policies(), (
            boost::python::arg("pixel"),
            boost::python::arg("panel") = 0,
            boost::python::arg("zstart") = 0,
            boost::python::arg("twod") = false)))
        .def("__init__", make_constructor(
          from_labels<2, FloatType>,
          default_call_policies(), (
            boost::python::arg("labels"),
            boost::python::arg("panel") = 0,
            boost::python::arg("zstart") = 0)))
        .def("__init__", make_constructor(
          from_labels<3, FloatType>,
          default_call_policies(), (
            boost::python::arg("labels"),
            boost::python::arg("panel") = 0)))
        .def("__init__", make_constructor(
          from_pixel_labeller<FloatType>,
          default_call_policies(), (
            boost::python::arg("labels"),
            boost::python::arg("panel") = 0)))
        .def("__init__", make_constructor(
          from_panel_and_bbox<FloatType>,
          default_call_policies(), (
            boost::python::arg("panel"),
            boost::python::arg("bbox"),
            boost::python::arg("allocate")=false,
            boost::python::arg("flatten")=false)))
        .def("allocate",
          &allocate<FloatType>)
        .def("allocate_with_value",
          &allocate_with_value<FloatType>)
        .def("deallocate",
          &deallocate<FloatType>)
        .def("is_consistent",
          &is_consistent<FloatType>)
        .def("panels",
          &panels<FloatType>)
        .def("bounding_boxes",
          &bounding_boxes<FloatType>)
        .def("count_mask_values",
          &count_mask_values<FloatType>)
        .def("is_bbox_within_image_volume",
          &is_bbox_within_image_volume<FloatType>, (
            boost::python::arg("image_size"),
            boost::python::arg("scan_range")))
        .def("does_bbox_contain_bad_pixels",
          &does_bbox_contain_bad_pixels<FloatType>, (
            boost::python::arg("mask")))
        .def("peak_coordinates",
          &peak_coordinates<FloatType>)
       .def("centroid_all",
          &centroid_all<FloatType>)
        .def("centroid_masked",
          &centroid_masked<FloatType>)
        .def("centroid_valid",
          &centroid_valid<FloatType>)
        .def("centroid_foreground",
          &centroid_foreground<FloatType>)
        .def("centroid_strong",
          &centroid_strong<FloatType>)
        .def("centroid_all_minus_background",
          &centroid_all_minus_background<FloatType>)
        .def("centroid_masked_minus_background",
          &centroid_masked_minus_background<FloatType>)
        .def("centroid_valid_minus_background",
          &centroid_valid_minus_background<FloatType>)
        .def("centroid_foreground_minus_background",
          &centroid_foreground_minus_background<FloatType>)
        .def("centroid_strong_minus_background",
          &centroid_strong_minus_background<FloatType>)
        .def("summed_intensity",
          &summed_intensity<FloatType>)
        .def("mean_background",
          &mean_background<FloatType>)
        .def("flatten", &flatten<FloatType>)
        .def_pickle(flex_pickle_double_buffered<shoebox_type,
          shoebox_to_string<FloatType>,
          shoebox_from_string<FloatType> >());
  }

  void export_flex_shoebox() {
    flex_shoebox_wrapper<ProfileFloatType>("shoebox");
  }

}}} // namespace dials::af::boost_python
