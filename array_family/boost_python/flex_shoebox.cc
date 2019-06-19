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
#include <cctbx/miller.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/pixel_list.h>
#include <dials/model/data/observation.h>
#include <dials/algorithms/image/connected_components/connected_components.h>
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using af::int2;
  using af::int6;
  using af::small;
  using dials::algorithms::LabelImageStack;
  using dials::algorithms::LabelPixels;
  using dials::algorithms::PixelToMillerIndex;
  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using dials::model::Centroid;
  using dials::model::Foreground;
  using dials::model::Intensity;
  using dials::model::Observation;
  using dials::model::PixelListLabeller;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dxtbx::model::BeamBase;
  using dxtbx::model::CrystalBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using scitbx::vec3;

  /**
   * Construct from an array of panels and bounding boxes.
   */
  template <typename FloatType>
  typename af::flex<Shoebox<FloatType> >::type *from_panel_and_bbox(
    const af::const_ref<std::size_t> panel,
    const af::const_ref<int6> bbox,
    bool allocate,
    bool flatten) {
    DIALS_ASSERT(panel.size() == bbox.size());
    af::shared<Shoebox<FloatType> > result(panel.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = Shoebox<FloatType>(panel[i], bbox[i], flatten);
      if (allocate) {
        result[i].allocate();
      }
    }
    return new typename af::flex<Shoebox<FloatType> >::type(
      result, af::flex_grid<>(result.size()));
  }

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <typename FloatType>
  class PixelListShoeboxCreator {
  public:
    PixelListShoeboxCreator(const PixelListLabeller &pixel,
                            std::size_t panel,
                            std::size_t zstart,
                            bool twod,
                            std::size_t min_pixels,
                            std::size_t max_pixels,
                            bool find_hot_pixels) {
      // Check the input
      DIALS_ASSERT(min_pixels > 0);
      DIALS_ASSERT(max_pixels > min_pixels);

      // Get the stuff from the label struct
      af::shared<int> labels = twod ? pixel.labels_2d() : pixel.labels_3d();
      af::shared<double> values = pixel.values();
      af::shared<vec3<int> > coords = pixel.coords();
      DIALS_ASSERT(labels.size() == values.size());
      DIALS_ASSERT(labels.size() == coords.size());

      // Get the number of labels and allocate the array
      std::size_t num = af::max(labels.const_ref()) + 1;
      result_ = af::shared<Shoebox<FloatType> >(num, Shoebox<FloatType>());
      spot_size_ = af::shared<std::size_t>(num, 0);

      // Initialise the bboxes
      int xsize = pixel.size()[1];
      int ysize = pixel.size()[0];
      int2 minmaxz = pixel.frame_range();
      for (std::size_t i = 0; i < result_.size(); ++i) {
        result_[i].panel = panel;
        result_[i].bbox[0] = xsize;
        result_[i].bbox[1] = 0;
        result_[i].bbox[2] = ysize;
        result_[i].bbox[3] = 0;
        result_[i].bbox[4] = minmaxz[1];
        result_[i].bbox[5] = minmaxz[0];
      }

      // Set the shoeboxes
      std::vector<std::size_t> num_pixels(num);
      for (std::size_t i = 0; i < labels.size(); ++i) {
        int l = labels[i];
        vec3<int> c = coords[i];
        DIALS_ASSERT(l < num_pixels.size());
        DIALS_ASSERT(c[2] < xsize && c[2] >= 0);
        DIALS_ASSERT(c[1] < ysize && c[1] >= 0);
        DIALS_ASSERT(c[0] < minmaxz[1] && c[0] >= minmaxz[0]);
        if (c[2] < result_[l].bbox[0]) result_[l].bbox[0] = c[2];
        if (c[2] >= result_[l].bbox[1]) result_[l].bbox[1] = c[2] + 1;
        if (c[1] < result_[l].bbox[2]) result_[l].bbox[2] = c[1];
        if (c[1] >= result_[l].bbox[3]) result_[l].bbox[3] = c[1] + 1;
        if (c[0] < result_[l].bbox[4]) result_[l].bbox[4] = c[0];
        if (c[0] >= result_[l].bbox[5]) result_[l].bbox[5] = c[0] + 1;
        num_pixels[l]++;
      }

      // Copy the spot size
      std::copy(num_pixels.begin(), num_pixels.end(), spot_size_.begin());

      // Allocate all the arrays
      for (std::size_t i = 0; i < result_.size(); ++i) {
        if (min_pixels <= num_pixels[i] && num_pixels[i] <= max_pixels) {
          result_[i].allocate();
        }
      }

      // Set all the mask and data points
      for (std::size_t i = 0; i < labels.size(); ++i) {
        int l = labels[i];
        if (result_[l].is_allocated()) {
          FloatType v = values[i];
          vec3<int> c = coords[i];
          int ii = c[2] - result_[l].bbox[0];
          int jj = c[1] - result_[l].bbox[2];
          int kk = c[0] - result_[l].bbox[4];
          DIALS_ASSERT(ii >= 0 && jj >= 0 && kk >= 0);
          DIALS_ASSERT(ii < result_[l].xsize());
          DIALS_ASSERT(jj < result_[l].ysize());
          DIALS_ASSERT(kk < result_[l].zsize());
          result_[l].data(kk, jj, ii) = v;
          result_[l].mask(kk, jj, ii) = Valid | Foreground;
        }
      }

      // Find hot pixels
      if (find_hot_pixels) {
        int first_frame = minmaxz[0];
        int last_frame = minmaxz[1];
        af::versa<int, af::c_grid<2> > hot_mask(af::c_grid<2>(ysize, xsize),
                                                first_frame - 1);
        for (std::size_t i = 0; i < coords.size(); ++i) {
          vec3<int> c = coords[i];
          if (c[0] == first_frame) {
            hot_mask(c[1], c[2]) = c[0];
          } else {
            if (hot_mask(c[1], c[2]) != c[0] - 1) {
              hot_mask(c[1], c[2]) = first_frame - 1;
            } else {
              hot_mask(c[1], c[2]) = c[0];
            }
          }
        }
        for (std::size_t i = 0; i < hot_mask.size(); ++i) {
          if (hot_mask[i] == last_frame - 1) {
            hot_pixels_.push_back(i);
          }
        }
      }

      // Shift bbox z start position
      for (std::size_t i = 0; i < result_.size(); ++i) {
        result_[i].bbox[4] += zstart;
        result_[i].bbox[5] += zstart;
      }
    }

    af::shared<Shoebox<FloatType> > result() const {
      DIALS_ASSERT(result_.size() == spot_size_.size());
      return result_;
    }

    af::shared<std::size_t> spot_size() const {
      DIALS_ASSERT(result_.size() == spot_size_.size());
      return spot_size_;
    }

    af::shared<std::size_t> hot_pixels() const {
      return hot_pixels_;
    }

  private:
    af::shared<Shoebox<FloatType> > result_;
    af::shared<std::size_t> spot_size_;
    af::shared<std::size_t> hot_pixels_;
  };

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <std::size_t DIM, typename FloatType>
  typename af::flex<Shoebox<FloatType> >::type *from_labels(
    const LabelImageStack<DIM> &label,
    std::size_t panel,
    std::size_t zstart) {
    // Get the stuff from the label struct
    af::shared<int> labels = label.labels();
    af::shared<int> values = label.values();
    af::shared<vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = af::max(labels.const_ref()) + 1;
    af::shared<Shoebox<FloatType> > result(num, Shoebox<FloatType>());

    // Initialise the bboxes
    int xsize = label.size()[1];
    int ysize = label.size()[0];
    int zsize = label.num_images();
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].bbox[0] = xsize;
      result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize;
      result[i].bbox[3] = 0;
      result[i].bbox[4] = zsize;
      result[i].bbox[5] = 0;
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] < result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] < result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] < result[l].bbox[4]) result[l].bbox[4] = c[0];
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
      result[l].data(kk, jj, ii) = v;
      result[l].mask(kk, jj, ii) = Valid | Foreground;
    }

    // Shift bbox z start position
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].bbox[4] += zstart;
      result[i].bbox[5] += zstart;
    }

    // Return the array
    return new
      typename af::flex<Shoebox<FloatType> >::type(result, af::flex_grid<>(num));
  }

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  template <typename FloatType>
  typename af::flex<Shoebox<FloatType> >::type *from_pixel_labeller(
    const LabelPixels &label,
    std::size_t panel) {
    // Get the stuff from the label struct
    af::shared<int> labels = label.labels();
    af::shared<int> values = label.values();
    af::shared<vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = af::max(labels.const_ref()) + 1;
    af::shared<Shoebox<FloatType> > result(num, Shoebox<FloatType>());

    // Initialise the bboxes
    int xsize = label.size()[2];
    int ysize = label.size()[1];
    int zsize = label.size()[0];
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].bbox[0] = xsize;
      result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize;
      result[i].bbox[3] = 0;
      result[i].bbox[4] = zsize;
      result[i].bbox[5] = 0;
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] < result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] < result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] < result[l].bbox[4]) result[l].bbox[4] = c[0];
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
      result[l].data(kk, jj, ii) = (double)v;
      result[l].mask(kk, jj, ii) = Valid | Foreground;
    }

    // Return the array
    return new
      typename af::flex<Shoebox<FloatType> >::type(result, af::flex_grid<>(num));
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
  shared<bool> is_consistent(const const_ref<Shoebox<FloatType> > &a) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_consistent();
    }
    return result;
  }

  /**
   * Check if the arrays are allocated
   */
  template <typename FloatType>
  shared<bool> is_allocated(const const_ref<Shoebox<FloatType> > &a) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_allocated();
    }
    return result;
  }

  /**
   * Check if the bounding box has points outside the image range.
   */
  template <typename FloatType>
  shared<bool> is_bbox_within_image_volume(const const_ref<Shoebox<FloatType> > &a,
                                           int2 image_size,
                                           int2 scan_range) {
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
  shared<bool> does_bbox_contain_bad_pixels(const const_ref<Shoebox<FloatType> > &a,
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
  shared<int> count_mask_values(const const_ref<Shoebox<FloatType> > &a, int code) {
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
  shared<vec3<double> > peak_coordinates(ref<Shoebox<FloatType> > a) {
    shared<vec3<double> > result(a.size(), af::init_functor_null<vec3<double> >());
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
  shared<int6> bounding_boxes(const const_ref<Shoebox<FloatType> > &a) {
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
  shared<std::size_t> panels(const const_ref<Shoebox<FloatType> > &a) {
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
  af::shared<Centroid> centroid_all(const const_ref<Shoebox<FloatType> > &a) {
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
  af::shared<Centroid> centroid_masked(const const_ref<Shoebox<FloatType> > &a,
                                       int code) {
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
  af::shared<Centroid> centroid_valid(const const_ref<Shoebox<FloatType> > &a) {
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
  af::shared<Centroid> centroid_foreground(const const_ref<Shoebox<FloatType> > &a) {
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
  af::shared<Centroid> centroid_strong(const const_ref<Shoebox<FloatType> > &a) {
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
    const const_ref<Shoebox<FloatType> > &a) {
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
    const const_ref<Shoebox<FloatType> > &a,
    int code) {
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
    const const_ref<Shoebox<FloatType> > &a) {
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
    const const_ref<Shoebox<FloatType> > &a) {
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
    const const_ref<Shoebox<FloatType> > &a) {
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
  af::shared<Intensity> summed_intensity(const const_ref<Shoebox<FloatType> > &a) {
    af::shared<Intensity> result(a.size(), Intensity());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].summed_intensity();
    }
    return result;
  }

  /**
   * Get a list of intensities
   */
  template <typename FloatType>
  af::shared<Intensity> bayesian_intensity(const const_ref<Shoebox<FloatType> > &a) {
    af::shared<Intensity> result(a.size(), Intensity());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = a[i].bayesian_intensity();
    }
    return result;
  }

  /**
   * Get the mean background.
   */
  template <typename FloatType>
  af::shared<double> mean_background(const const_ref<Shoebox<FloatType> > &a) {
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
   * Get the mean background.
   */
  template <typename FloatType>
  af::shared<double> mean_modelled_background(const const_ref<Shoebox<FloatType> > &a) {
    af::shared<double> result(a.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      af::versa<FloatType, af::c_grid<3> > data = a[i].background;
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
  void flatten(ref<Shoebox<FloatType> > self) {
    for (std::size_t i = 0; i < self.size(); ++i) {
      self[i].flatten();
    }
  }

  /**
   * Apply the shoebox mask to the background mask
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > apply_background_mask(
    const af::const_ref<Shoebox<FloatType> > &self,
    int frame,
    std::size_t num_panels,
    int2 size) {
    DIALS_ASSERT(num_panels > 0);
    DIALS_ASSERT(size[0] > 0);
    DIALS_ASSERT(size[1] > 0);
    DIALS_ASSERT(num_panels == 1);
    af::versa<bool, af::c_grid<2> > mask(af::c_grid<2>(size[0], size[1]), true);
    int height = size[0];
    int width = size[1];
    for (std::size_t s = 0; s < self.size(); ++s) {
      // Get stuff from the shoebox
      std::size_t p = self[s].panel;
      int x0 = self[s].bbox[0];
      int x1 = self[s].bbox[1];
      int y0 = self[s].bbox[2];
      int y1 = self[s].bbox[3];
      int z0 = self[s].bbox[4];
      DIALS_ASSERT(p == 0);

      // Get the shoebox mask
      af::const_ref<int, af::c_grid<3> > sbox_mask = self[s].mask.const_ref();

      // Make sure bbox range is ok
      int x00 = x0 >= 0 ? x0 : 0;
      int y00 = y0 >= 0 ? y0 : 0;
      int x11 = x1 <= width ? x1 : width;
      int y11 = y1 <= height ? y1 : height;

      // Set the mask
      int k = frame - z0;
      DIALS_ASSERT(k >= 0);
      DIALS_ASSERT(k < sbox_mask.accessor()[0]);
      for (std::size_t y = y00; y < y11; ++y) {
        for (std::size_t x = x00; x < x11; ++x) {
          std::size_t j = y - y0;
          std::size_t i = x - x0;
          DIALS_ASSERT(j < sbox_mask.accessor()[1]);
          DIALS_ASSERT(i < sbox_mask.accessor()[2]);
          mask(y, x) &= ((sbox_mask(k, j, i) & Background) != 0);
        }
      }
    }
    return mask;
  }

  /**
   * Apply the data, mask and background to the shoebox
   */
  template <typename FloatType>
  void apply_pixel_data(af::shared<Shoebox<FloatType> > self,
                        const af::const_ref<double, af::c_grid<2> > &data,
                        const af::const_ref<double, af::c_grid<2> > &background,
                        const af::const_ref<bool, af::c_grid<2> > &mask,
                        int frame,
                        std::size_t num_panels) {
    DIALS_ASSERT(num_panels > 0);
    DIALS_ASSERT(num_panels == 1);
    int height = background.accessor()[0];
    int width = background.accessor()[1];
    for (std::size_t s = 0; s < self.size(); ++s) {
      // Check shoebox
      DIALS_ASSERT(self[s].is_consistent());

      // Get stuff from the shoebox
      std::size_t p = self[s].panel;
      int x0 = self[s].bbox[0];
      int x1 = self[s].bbox[1];
      int y0 = self[s].bbox[2];
      int y1 = self[s].bbox[3];
      int z0 = self[s].bbox[4];
      DIALS_ASSERT(p == 0);

      // Get the shoebox mask
      af::ref<FloatType, af::c_grid<3> > sbox_data = self[s].data.ref();
      af::ref<FloatType, af::c_grid<3> > sbox_bgrd = self[s].background.ref();
      af::ref<int, af::c_grid<3> > sbox_mask = self[s].mask.ref();

      // Make sure bbox range is ok
      int x00 = x0 >= 0 ? x0 : 0;
      int y00 = y0 >= 0 ? y0 : 0;
      int x11 = x1 <= width ? x1 : width;
      int y11 = y1 <= height ? y1 : height;

      // Set the mask
      int k = frame - z0;
      DIALS_ASSERT(k >= 0);
      DIALS_ASSERT(k < sbox_data.accessor()[0]);
      for (std::size_t y = y00; y < y11; ++y) {
        for (std::size_t x = x00; x < x11; ++x) {
          std::size_t j = y - y0;
          std::size_t i = x - x0;
          DIALS_ASSERT(j < sbox_data.accessor()[1]);
          DIALS_ASSERT(i < sbox_data.accessor()[2]);
          sbox_data(k, j, i) = data(y, x);
          sbox_bgrd(k, j, i) = background(y, x);
          if (mask(y, x) == false) {
            sbox_mask(k, j, i) &= ~Valid;
          } else {
            sbox_mask(k, j, i) |= Valid;
          }
        }
      }
    }
  }

  /**
   * Do a pixel labelling filter to set any pixel closer to another pixel to
   * background
   */
  template <typename FloatType>
  bool mask_neighbouring_single(Shoebox<FloatType> &self,
                                cctbx::miller::index<> hkl,
                                const PixelToMillerIndex &compute_miller_index) {
    bool modified = false;
    int mask_code = Valid | Background;
    int x0 = self.bbox[0];
    int x1 = self.bbox[1];
    int y0 = self.bbox[2];
    int y1 = self.bbox[3];
    int z0 = self.bbox[4];
    int z1 = self.bbox[5];
    DIALS_ASSERT(x0 < x1);
    DIALS_ASSERT(y0 < y1);
    DIALS_ASSERT(z0 < z1);
    DIALS_ASSERT(self.is_consistent());
    std::size_t xsize = self.xsize();
    std::size_t ysize = self.ysize();
    std::size_t zsize = self.zsize();
    for (std::size_t z = 0; z < zsize; ++z) {
      for (std::size_t y = 0; y < ysize; ++y) {
        for (std::size_t x = 0; x < xsize; ++x) {
          double z1 = z0 + z;
          double y1 = y0 + y + 0.5;
          double x1 = x0 + x + 0.5;
          vec3<double> pixel_hkl1 = compute_miller_index.h(self.panel, x1, y1, z1);
          vec3<double> pixel_hkl2 = compute_miller_index.h(self.panel, x1, y1, z1 + 1);
          int h1 = (int)std::floor(pixel_hkl1[0] + 0.5);
          int k1 = (int)std::floor(pixel_hkl1[1] + 0.5);
          int l1 = (int)std::floor(pixel_hkl1[2] + 0.5);
          int h2 = (int)std::floor(pixel_hkl2[0] + 0.5);
          int k2 = (int)std::floor(pixel_hkl2[1] + 0.5);
          int l2 = (int)std::floor(pixel_hkl2[2] + 0.5);
          if (h1 != hkl[0] || k1 != hkl[1] || l1 != hkl[2] || h2 != hkl[0]
              || k2 != hkl[1] || l2 != hkl[2]) {
            self.mask(z, y, x) = mask_code;
            modified = true;
          }
        }
      }
    }
    return modified;
  }

  /**
   * Do a pixel labelling filter to set any pixel closer to another pixel to
   * background
   */
  template <typename FloatType>
  af::shared<bool> mask_neighbouring(af::ref<Shoebox<FloatType> > self,
                                     af::const_ref<cctbx::miller::index<> > hkl,
                                     const BeamBase &beam,
                                     const Detector &detector,
                                     const Goniometer &goniometer,
                                     const Scan &scan,
                                     const CrystalBase &crystal) {
    DIALS_ASSERT(self.size() == hkl.size());
    af::shared<bool> modified(self.size());
    PixelToMillerIndex compute_miller_index(beam, detector, goniometer, scan, crystal);
    for (std::size_t i = 0; i < self.size(); ++i) {
      modified[i] = mask_neighbouring_single(self[i], hkl[i], compute_miller_index);
    }
    return modified;
  }

  /**
   * A class to convert the shoebox class to a string for pickling
   */
  template <typename FloatType>
  struct shoebox_to_string : pickle_double_buffered::to_string {
    using pickle_double_buffered::to_string::operator<<;

    typedef Shoebox<FloatType> shoebox_type;

    /** Initialise with the version for checking */
    shoebox_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    shoebox_to_string &operator<<(const shoebox_type &val) {
      *this << val.panel << val.bbox[0] << val.bbox[1] << val.bbox[2] << val.bbox[3]
            << val.bbox[4] << val.bbox[5];

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
  struct shoebox_from_string : pickle_double_buffered::from_string {
    using pickle_double_buffered::from_string::operator>>;

    typedef Shoebox<FloatType> shoebox_type;

    /** Initialise the class with the string. Get the version and check */
    shoebox_from_string(const char *str_ptr)
        : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    shoebox_from_string &operator>>(shoebox_type &val) {
      *this >> val.panel >> val.bbox[0] >> val.bbox[1] >> val.bbox[2] >> val.bbox[3]
        >> val.bbox[4] >> val.bbox[5];

      val.data = profile_from_string<versa<FloatType, c_grid<3> > >();
      val.mask = profile_from_string<versa<int, c_grid<3> > >();
      val.background = profile_from_string<versa<FloatType, c_grid<3> > >();

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
  typename scitbx::af::boost_python::
    flex_wrapper<Shoebox<FloatType>, return_internal_reference<> >::class_f_t
    flex_shoebox_wrapper(const char *name) {
    typedef Shoebox<FloatType> shoebox_type;

    return scitbx::af::boost_python::
      flex_wrapper<shoebox_type, return_internal_reference<> >::plain(name)
        /* .def("__init__", make_constructor( */
        /*   from_pixel_list<FloatType>, */
        /*   default_call_policies(), ( */
        /*     boost::python::arg("pixel"), */
        /*     boost::python::arg("panel") = 0, */
        /*     boost::python::arg("zstart") = 0, */
        /*     boost::python::arg("twod") = false, */
        /*     boost::python::arg("min_pixels") = 1, */
        /*     boost::python::arg("max_pixels") = 20))) */
        .def("__init__",
             make_constructor(from_labels<2, FloatType>,
                              default_call_policies(),
                              (boost::python::arg("labels"),
                               boost::python::arg("panel") = 0,
                               boost::python::arg("zstart") = 0)))
        .def("__init__",
             make_constructor(
               from_labels<3, FloatType>,
               default_call_policies(),
               (boost::python::arg("labels"), boost::python::arg("panel") = 0)))
        .def("__init__",
             make_constructor(
               from_pixel_labeller<FloatType>,
               default_call_policies(),
               (boost::python::arg("labels"), boost::python::arg("panel") = 0)))
        .def("__init__",
             make_constructor(from_panel_and_bbox<FloatType>,
                              default_call_policies(),
                              (boost::python::arg("panel"),
                               boost::python::arg("bbox"),
                               boost::python::arg("allocate") = false,
                               boost::python::arg("flatten") = false)))
        .def("allocate", &allocate<FloatType>)
        .def("allocate_with_value", &allocate_with_value<FloatType>)
        .def("deallocate", &deallocate<FloatType>)
        .def("is_consistent", &is_consistent<FloatType>)
        .def("is_allocated", &is_allocated<FloatType>)
        .def("panels", &panels<FloatType>)
        .def("bounding_boxes", &bounding_boxes<FloatType>)
        .def("count_mask_values", &count_mask_values<FloatType>)
        .def("is_bbox_within_image_volume",
             &is_bbox_within_image_volume<FloatType>,
             (boost::python::arg("image_size"), boost::python::arg("scan_range")))
        .def("does_bbox_contain_bad_pixels",
             &does_bbox_contain_bad_pixels<FloatType>,
             (boost::python::arg("mask")))
        .def("peak_coordinates", &peak_coordinates<FloatType>)
        .def("centroid_all", &centroid_all<FloatType>)
        .def("centroid_masked", &centroid_masked<FloatType>)
        .def("centroid_valid", &centroid_valid<FloatType>)
        .def("centroid_foreground", &centroid_foreground<FloatType>)
        .def("centroid_strong", &centroid_strong<FloatType>)
        .def("centroid_all_minus_background", &centroid_all_minus_background<FloatType>)
        .def("centroid_masked_minus_background",
             &centroid_masked_minus_background<FloatType>)
        .def("centroid_valid_minus_background",
             &centroid_valid_minus_background<FloatType>)
        .def("centroid_foreground_minus_background",
             &centroid_foreground_minus_background<FloatType>)
        .def("centroid_strong_minus_background",
             &centroid_strong_minus_background<FloatType>)
        .def("bayesian_intensity", &bayesian_intensity<FloatType>)
        .def("summed_intensity", &summed_intensity<FloatType>)
        .def("mean_background", &mean_background<FloatType>)
        .def("mean_modelled_background", &mean_modelled_background<FloatType>)
        .def("flatten", &flatten<FloatType>)
        .def("apply_background_mask", &apply_background_mask<FloatType>)
        .def("apply_pixel_data", &apply_pixel_data<FloatType>)
        .def("mask_neighbouring", &mask_neighbouring<FloatType>)
        .def_pickle(flex_pickle_double_buffered<shoebox_type,
                                                shoebox_to_string<FloatType>,
                                                shoebox_from_string<FloatType> >());
  }

  void export_flex_shoebox() {
    flex_shoebox_wrapper<ProfileFloatType>("shoebox");

    class_<PixelListShoeboxCreator<ProfileFloatType> >("PixelListShoeboxCreator",
                                                       no_init)
      .def(init<const PixelListLabeller &,
                std::size_t,
                std::size_t,
                bool,
                std::size_t,
                std::size_t,
                bool>((boost::python::arg("pixel"),
                       boost::python::arg("panel") = 0,
                       boost::python::arg("zstart") = 0,
                       boost::python::arg("twod") = false,
                       boost::python::arg("min_pixels") = 1,
                       boost::python::arg("max_pixels") = 20,
                       boost::python::arg("find_hot_pixels") = false)))
      .def("result", &PixelListShoeboxCreator<ProfileFloatType>::result)
      .def("spot_size", &PixelListShoeboxCreator<ProfileFloatType>::spot_size)
      .def("hot_pixels", &PixelListShoeboxCreator<ProfileFloatType>::hot_pixels);
  }

}}}  // namespace dials::af::boost_python
