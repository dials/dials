/*
 * processor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROCESSOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PROCESSOR_H

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <ctime>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dxtbx/array_family/flex_table_suite.h>
#include <dials/array_family/boost_python/reflection_table_suite.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <chrono>

namespace dials { namespace algorithms {

  using model::Image;
  using model::Shoebox;
  using model::Valid;

  /**
   * The cctbx build system is too messed up to figure out how to build
   * boost::system need by boost::chrono. Therefore just use system clock
   * to get a timestamp in ms
   */
  double timestamp() {
    return ((double)clock()) / ((double)CLOCKS_PER_SEC);
  }

  /**
   * A base class for executor callbacks
   */
  class Executor {
  public:
    virtual void process(int, af::reflection_table) = 0;
  };

  /**
   * A class to extract shoebox pixels from images
   */
  class ShoeboxProcessor {
  public:
    /**
     * Initialise the index array. Determine which reflections are recorded on
     * each frame and panel ahead of time to enable quick lookup of the
     * reflections to be written to when processing each image.
     */
    ShoeboxProcessor(af::reflection_table data,
                     std::size_t npanels,
                     int frame0,
                     int frame1,
                     bool save)
        : data_(data),
          extract_time_(0.0),
          process_time_(0.0),
          save_(save),
          npanels_(npanels),
          frame0_(frame0),
          frame1_(frame1),
          frame_(frame0),
          nframes_(frame1 - frame0) {
      DIALS_ASSERT(frame0_ < frame1_);
      DIALS_ASSERT(npanels_ > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.size() > 0);
      af::const_ref<Shoebox<>> shoebox = data["shoebox"];
      std::size_t size = nframes_ * npanels_;
      std::vector<std::size_t> num(size, 0);
      std::vector<std::size_t> count(size, 0);
      flatten_ = shoebox[0].flat;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        DIALS_ASSERT(shoebox[i].flat == flatten_);
        DIALS_ASSERT(shoebox[i].is_allocated() == false);
        DIALS_ASSERT(shoebox[i].bbox[1] > shoebox[i].bbox[0]);
        DIALS_ASSERT(shoebox[i].bbox[3] > shoebox[i].bbox[2]);
        DIALS_ASSERT(shoebox[i].bbox[5] > shoebox[i].bbox[4]);
        for (int z = shoebox[i].bbox[4]; z < shoebox[i].bbox[5]; ++z) {
          std::size_t j = shoebox[i].panel + (z - frame0_) * npanels_;
          DIALS_ASSERT(j < num.size());
          num[j]++;
        }
      }
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));
      indices_.resize(offset_.back());
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        for (int z = shoebox[i].bbox[4]; z < shoebox[i].bbox[5]; ++z) {
          std::size_t j = shoebox[i].panel + (z - frame0_) * npanels_;
          std::size_t k = offset_[j] + count[j];
          DIALS_ASSERT(j < count.size());
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = i;
          count[j]++;
        }
      }
      DIALS_ASSERT(count == num);
    }

    /**
     * Extract the pixels from the image and copy to the relevant shoeboxes.
     * @param image The image to process
     * @param frame The current image frame
     */
    template <typename T>
    void next(const Image<T>& image, Executor& executor) {
      using dials::af::boost_python::reflection_table_suite::select_rows_index;
      using dxtbx::af::flex_table_suite::set_selected_rows_index;
      typedef Shoebox<>::float_type float_type;
      typedef af::ref<float_type, af::c_grid<3>> sbox_data_type;
      typedef af::ref<int, af::c_grid<3>> sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);

      // Get the initial time
      double start_time = timestamp();

      // For each image, extract shoeboxes of reflections recorded.
      // Allocate data where necessary
      af::ref<Shoebox<>> shoebox = data_["shoebox"];
      af::shared<std::size_t> process_indices;
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2>> data = image.data(p);
        af::const_ref<bool, af::c_grid<2>> mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox.size());
          Shoebox<>& sbox = shoebox[ind[i]];
          if (frame_ == sbox.bbox[4]) {
            DIALS_ASSERT(sbox.is_allocated() == false);
            sbox.allocate();
          }
          int6 b = sbox.bbox;
          sbox_data_type sdata = sbox.data.ref();
          sbox_mask_type smask = sbox.mask.ref();
          DIALS_ASSERT(b[1] > b[0]);
          DIALS_ASSERT(b[3] > b[2]);
          DIALS_ASSERT(b[5] > b[4]);
          DIALS_ASSERT(frame_ >= b[4] && frame_ < b[5]);
          int x0 = b[0];
          int x1 = b[1];
          int y0 = b[2];
          int y1 = b[3];
          int z0 = b[4];
          int xs = x1 - x0;
          int ys = y1 - y0;
          int z = frame_ - z0;
          int yi = (int)data.accessor()[0];
          int xi = (int)data.accessor()[1];
          int xb = x0 >= 0 ? 0 : std::abs(x0);
          int yb = y0 >= 0 ? 0 : std::abs(y0);
          int xe = x1 <= xi ? xs : xs - (x1 - xi);
          int ye = y1 <= yi ? ys : ys - (y1 - yi);
          if (yb >= ye || xb >= xe) {
            continue;
          }
          DIALS_ASSERT(yb >= 0 && ye <= ys);
          DIALS_ASSERT(xb >= 0 && xe <= xs);
          DIALS_ASSERT(yb + y0 >= 0 && ye + y0 <= yi);
          DIALS_ASSERT(xb + x0 >= 0 && xe + x0 <= xi);
          DIALS_ASSERT(sbox.is_consistent());
          if (flatten_) {
            for (std::size_t y = yb; y < ye; ++y) {
              for (std::size_t x = xb; x < xe; ++x) {
                sdata(0, y, x) += data(y + y0, x + x0);
                bool sv = smask(0, y, x) & Valid;
                bool mv = mask(y + y0, x + x0);
                smask(0, y, x) = (mv && (z == 0 ? true : sv) ? Valid : 0);
              }
            }
          } else {
            for (std::size_t y = yb; y < ye; ++y) {
              for (std::size_t x = xb; x < xe; ++x) {
                sdata(z, y, x) = data(y + y0, x + x0);
                smask(z, y, x) = mask(y + y0, x + x0) ? Valid : 0;
              }
            }
          }
          if (frame_ == sbox.bbox[5] - 1) {
            process_indices.push_back(ind[i]);
          }
        }
      }

      // Update timing info
      double end_time = timestamp();
      extract_time_ += end_time - start_time;

      // Process all the reflections and set the reflections
      if (process_indices.size() > 0) {
        double start_time = timestamp();
        af::const_ref<std::size_t> ind = process_indices.const_ref();
        af::reflection_table reflections = select_rows_index(data_, ind);
        executor.process(frame_, reflections);
        set_selected_rows_index(data_, ind, reflections);
        if (!save_) {
          for (std::size_t i = 0; i < ind.size(); ++i) {
            shoebox[ind[i]].deallocate();
          }
        }
        double end_time = timestamp();
        process_time_ += end_time - start_time;
      }

      // Update the frame counter
      frame_++;
    }

    /**
     * Extract the pixels from the image and copy to the relevant shoeboxes.
     * @param image The image to process
     * @param frame The current image frame
     */
    template <typename T>
    void next_data_only(const Image<T>& image) {
      using dials::af::boost_python::reflection_table_suite::select_rows_index;
      using dxtbx::af::flex_table_suite::set_selected_rows_index;
      typedef Shoebox<>::float_type float_type;
      typedef af::ref<float_type, af::c_grid<3>> sbox_data_type;
      typedef af::ref<int, af::c_grid<3>> sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);

      // Get the initial time
      double start_time = timestamp();

      // For each image, extract shoeboxes of reflections recorded.
      // Allocate data where necessary
      af::ref<Shoebox<>> shoebox = data_["shoebox"];
      af::shared<std::size_t> process_indices;
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2>> data = image.data(p);
        af::const_ref<bool, af::c_grid<2>> mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox.size());
          Shoebox<>& sbox = shoebox[ind[i]];
          if (frame_ == sbox.bbox[4]) {
            DIALS_ASSERT(sbox.is_allocated() == false);
            sbox.allocate();
          }
          int6 b = sbox.bbox;
          sbox_data_type sdata = sbox.data.ref();
          sbox_mask_type smask = sbox.mask.ref();
          DIALS_ASSERT(b[1] > b[0]);
          DIALS_ASSERT(b[3] > b[2]);
          DIALS_ASSERT(b[5] > b[4]);
          DIALS_ASSERT(frame_ >= b[4] && frame_ < b[5]);
          int x0 = b[0];
          int x1 = b[1];
          int y0 = b[2];
          int y1 = b[3];
          int z0 = b[4];
          int xs = x1 - x0;
          int ys = y1 - y0;
          int z = frame_ - z0;
          int yi = (int)data.accessor()[0];
          int xi = (int)data.accessor()[1];
          int xb = x0 >= 0 ? 0 : std::abs(x0);
          int yb = y0 >= 0 ? 0 : std::abs(y0);
          int xe = x1 <= xi ? xs : xs - (x1 - xi);
          int ye = y1 <= yi ? ys : ys - (y1 - yi);
          if (yb >= ye || xb >= xe) {
            continue;
          }
          DIALS_ASSERT(yb >= 0 && ye <= ys);
          DIALS_ASSERT(xb >= 0 && xe <= xs);
          DIALS_ASSERT(yb + y0 >= 0 && ye + y0 <= yi);
          DIALS_ASSERT(xb + x0 >= 0 && xe + x0 <= xi);
          DIALS_ASSERT(sbox.is_consistent());
          if (flatten_) {
            for (std::size_t y = yb; y < ye; ++y) {
              for (std::size_t x = xb; x < xe; ++x) {
                sdata(0, y, x) += data(y + y0, x + x0);
                bool sv = smask(0, y, x) & Valid;
                bool mv = mask(y + y0, x + x0);
                smask(0, y, x) = (mv && (z == 0 ? true : sv) ? Valid : 0);
              }
            }
          } else {
            for (std::size_t y = yb; y < ye; ++y) {
              for (std::size_t x = xb; x < xe; ++x) {
                sdata(z, y, x) = data(y + y0, x + x0);
                smask(z, y, x) = mask(y + y0, x + x0) ? Valid : 0;
              }
            }
          }
          if (frame_ == sbox.bbox[5] - 1) {
            process_indices.push_back(ind[i]);
          }
        }
      }

      // Update timing info
      double end_time = timestamp();
      extract_time_ += end_time - start_time;

      // Update the frame counter
      frame_++;
    }
    /** @returns The first frame.  */
    int frame0() const {
      return frame0_;
    }

    /** @returns The last frame */
    int frame1() const {
      return frame1_;
    }

    /** @returns The current frame. */
    int frame() const {
      return frame_;
    }

    /** @returns The number of frames  */
    std::size_t nframes() const {
      return nframes_;
    }

    /** @returns The number of panels */
    std::size_t npanels() const {
      return npanels_;
    }

    /**
     * @returns Is the extraction finished.
     */
    bool finished() const {
      return frame_ == frame1_;
    }

    /**
     * @returns The extract time
     */
    double extract_time() const {
      return extract_time_;
    }

    /**
     * @returns The process time
     */
    double process_time() const {
      return process_time_;
    }

  private:
    /**
     * Get an index array specifying which reflections are recorded on a given
     * frame and panel.
     * @param frame The frame number
     * @param panel The panel number
     * @returns An array of indices
     */
    af::const_ref<std::size_t> indices(int frame, std::size_t panel) const {
      std::size_t j0 = panel + (frame - frame0_) * npanels_;
      DIALS_ASSERT(offset_.size() > 0);
      DIALS_ASSERT(j0 < offset_.size() - 1);
      std::size_t i0 = offset_[j0];
      std::size_t i1 = offset_[j0 + 1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - off;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t>(&indices_[off], num);
    }

    af::reflection_table data_;
    double extract_time_;
    double process_time_;
    bool flatten_;
    bool save_;
    std::size_t npanels_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
  };

  class ShoeboxProcessorV2 {
  public:
    ShoeboxProcessorV2(af::reflection_table data,
                       std::size_t npanels,
                       int frame0,
                       int frame1,
                       bool save,
                       const dxtbx::model::Scan& scan,
                       const dxtbx::model::BeamBase& beam,
                       const dxtbx::model::Goniometer& gonio,
                       const dxtbx::model::Detector& detector,
                       const double delta_b,
                       const double delta_m)
        : data_(data),
          extract_time_(0.0),
          process_time_(0.0),
          save_(save),
          npanels_(npanels),
          frame0_(frame0),
          frame1_(frame1),
          frame_(frame0),
          nframes_(frame1 - frame0),
          phi0_(scan.get_oscillation()[0]),
          dphi_(scan.get_oscillation()[1]),
          s0_(beam.get_s0()),
          m2_(gonio.get_rotation_axis()),
          detector_(detector),
          index0_(scan.get_array_range()[0]),
          index1_(scan.get_array_range()[1]) {
      delta_b_r2 = 1.0 / (std::pow(delta_b, 2));
      delta_m_r2 = 1.0 / (std::pow(delta_m, 2));
      DIALS_ASSERT(frame0_ < frame1_);
      DIALS_ASSERT(npanels_ > 0);
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.size() > 0);

      // Precalculate a few things.
      af::ref<Shoebox<>> shoebox = data_["shoebox"];
      af::ref<vec3<double>> s1_vec = data_["s1"];
      af::ref<vec3<double>> xyzcal_px = data_["xyzcal.px"];
      for (std::size_t p = 0; p < npanels_; ++p) {
        const dxtbx::model::Panel& panel = detector_[p];
        std::vector<profile_model::gaussian_rs::CoordinateSystem> csi;
        std::vector<double> atten_lengths_p;
        csi.reserve(data.size());
        atten_lengths_p.reserve(data.size());
        for (int i = 0; i < data_.size(); ++i) {
          vec3<double> s1 = s1_vec[i];
          vec3<double> xyzcal = xyzcal_px[i];
          double phi = phi0_ + (xyzcal[2] - index0_) * dphi_;
          profile_model::gaussian_rs::CoordinateSystem cs(m2_, s0_, s1, phi);
          csi.push_back(cs);
          vec2<double> shoebox_centroid_px = panel.get_ray_intersection_px(s1);
          double attenuation_length = panel.attenuation_length(shoebox_centroid_px);
          atten_lengths_p.push_back(attenuation_length);
        }
        coordinate_systems_[p] = csi;
        attenuation_lengths_[p] = atten_lengths_p;
      }

      std::size_t size = nframes_ * npanels_;
      std::vector<std::size_t> num(size, 0);
      std::vector<std::size_t> count(size, 0);
      flatten_ = shoebox[0].flat;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        DIALS_ASSERT(shoebox[i].flat == flatten_);
        DIALS_ASSERT(shoebox[i].bbox[1] > shoebox[i].bbox[0]);
        DIALS_ASSERT(shoebox[i].bbox[3] > shoebox[i].bbox[2]);
        DIALS_ASSERT(shoebox[i].bbox[5] > shoebox[i].bbox[4]);
        for (int z = shoebox[i].bbox[4]; z < shoebox[i].bbox[5]; ++z) {
          std::size_t j = shoebox[i].panel + (z - frame0_) * npanels_;
          DIALS_ASSERT(j < num.size());
          num[j]++;
        }
      }
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));
      indices_.resize(offset_.back());
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        for (int z = shoebox[i].bbox[4]; z < shoebox[i].bbox[5]; ++z) {
          std::size_t j = shoebox[i].panel + (z - frame0_) * npanels_;
          std::size_t k = offset_[j] + count[j];
          DIALS_ASSERT(j < count.size());
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = i;
          count[j]++;
        }
      }
      DIALS_ASSERT(count == num);
    }

    template <typename T>
    void next_data_only(const Image<T>& image) {
      using dials::af::boost_python::reflection_table_suite::select_rows_index;
      using dxtbx::af::flex_table_suite::set_selected_rows_index;
      typedef Shoebox<>::float_type float_type;
      typedef af::ref<float_type, af::c_grid<3>> sbox_data_type;
      typedef af::ref<int, af::c_grid<3>> sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);

      // Get the initial time
      double start_time = timestamp();

      // For each image, extract shoeboxes of reflections recorded.
      // Allocate data where necessary
      af::ref<Shoebox<>> shoebox = data_["shoebox"];
      double s0_length = s0_.length();

      for (std::size_t p = 0; p < image.npanels(); ++p) {
        std::vector<double> attenuation_lengths = attenuation_lengths_[p];
        const dxtbx::model::Panel& panel = detector_[p];
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2>> data = image.data(p);
        af::const_ref<bool, af::c_grid<2>> mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox.size());
          Shoebox<>& sbox = shoebox[ind[i]];
          int6 b = sbox.bbox;
          DIALS_ASSERT(b[1] > b[0]);
          DIALS_ASSERT(b[3] > b[2]);
          DIALS_ASSERT(b[5] > b[4]);
          DIALS_ASSERT(frame_ >= b[4] && frame_ < b[5]);
          int x0 = b[0];
          int x1 = b[1];
          int y0 = b[2];
          int y1 = b[3];
          int z0 = b[4];
          int xs = x1 - x0;
          int ys = y1 - y0;
          int z = frame_ - z0;
          int yi = (int)data.accessor()[0];
          int xi = (int)data.accessor()[1];
          // These bits make sure we only read within the image
          int xb = x0 >= 0 ? 0 : std::abs(x0);
          int yb = y0 >= 0 ? 0 : std::abs(y0);
          int xe = x1 <= xi ? xs : xs - (x1 - xi);
          int ye = y1 <= yi ? ys : ys - (y1 - yi);
          if (yb >= ye || xb >= xe) {
            continue;
          }
          DIALS_ASSERT(yb >= 0 && ye <= ys);
          DIALS_ASSERT(xb >= 0 && xe <= xs);
          DIALS_ASSERT(yb + y0 >= 0 && ye + y0 <= yi);
          DIALS_ASSERT(xb + x0 >= 0 && xe + x0 <= xi);

          double attenuation_length = attenuation_lengths[ind[i]];
          profile_model::gaussian_rs::CoordinateSystem cs =
            coordinate_systems_[p][ind[i]];

          // This implementation is slightly less efficient for this part, as it is
          // recalculated for each slice even though it is the same (this is not a
          // problem in the other integrators as it is done in the mask calculator). But
          // I think this calculation should be done in Kabsch coordinates, which would
          // be different for each slice anyway, so let's not worry about it.

          af::versa<double, af::c_grid<2>> dxy_array(af::c_grid<2>(ys + 1, xs + 1));
          for (int j2 = 0; j2 <= ys; ++j2) {
            for (int i2 = 0; i2 <= xs; ++i2) {
              vec2<double> gxy = cs.from_beam_vector(
                panel
                  .get_pixel_lab_coord(vec2<double>(x0 + i2, y0 + j2),
                                       attenuation_length)
                  .normalize()
                * s0_length);
              dxy_array(j2, i2) = (gxy[0] * gxy[0] + gxy[1] * gxy[1]) * delta_b_r2;
            }
          }

          // As we are combining data loading and mask calculation, we need to check
          // if parts of the bbox are outside of the image. This could be saved per
          // shoebox but is a pretty cheap calculation.
          bool all_in_image_bounds = true;
          af::versa<bool, af::c_grid<2>> in_image_array(af::c_grid<2>(ys, xs), true);
          if ((x0 < 0) || (y0 < 0) || ((x1 >= xi) || (y1 >= yi))) {
            all_in_image_bounds = false;
            for (int j3 = 0; j3 < ys; ++j3) {
              for (int i3 = 0; i3 < xs; ++i3) {
                if ((j3 + y0 < 0) || (i3 + x0 < 0) || (x0 + i3 >= xi)
                    || (y0 + j3 >= yi)) {
                  in_image_array(j3, i3) = false;
                }
              }
            }
          }

          // now loop through the shoebox pixels, test if they are in the image
          // and see if they are foreground, background, valid etc and
          // add them to the right quantity in the shoebox.
          for (int j3 = 0; j3 < ys; ++j3) {
            for (int i3 = 0; i3 < xs; ++i3) {
              bool this_in_image_bounds = all_in_image_bounds;
              if (!this_in_image_bounds) {  // if not all in bounds, check the
                                            // specific pixel
                this_in_image_bounds = in_image_array(j3, i3);
              }
              double dxy1 = dxy_array(j3, i3);
              double dxy2 = dxy_array(j3 + 1, i3);
              double dxy3 = dxy_array(j3, i3 + 1);
              double dxy4 = dxy_array(j3 + 1, i3 + 1);
              double dxy = std::min(std::min(dxy1, dxy2), std::min(dxy3, dxy4));

              if (dxy <= 1.0) {  // Is Foreground
                if (this_in_image_bounds) {
                  if (mask(j3 + y0, i3 + x0)) {  // Is Valid
                    sbox.total_intensity += data(j3 + y0, i3 + x0);
                    sbox.n_valid_fg += 1;
                  } else {
                    sbox.masked_image_pixel = true;
                    sbox.n_invalid_fg += 1;
                  }
                } else {
                  sbox.masked_image_pixel = true;
                  sbox.n_invalid_fg += 1;
                }
              } else {  // Is Background
                if (this_in_image_bounds) {
                  if (mask(j3 + y0, i3 + x0)) {  // Is Valid
                    sbox.n_valid_bg += 1;
                    int this_pixel = data(j3 + y0, i3 + x0);
                    if (auto search = sbox.background_hist.find(this_pixel);
                        search != sbox.background_hist.end()) {
                      sbox.background_hist[this_pixel] += 1;
                    } else {
                      sbox.background_hist[this_pixel] = 1;
                    }
                  } else {
                    sbox.n_invalid_bg += 1;
                  }
                } else {
                  sbox.n_invalid_bg += 1;
                }
              }
            }
          }
        }
      }

      // Update timing info
      double end_time = timestamp();
      extract_time_ += end_time - start_time;

      // Update the frame counter
      frame_++;
    }

    template <typename T>
    af::shared<double> finalise(af::reflection_table& data) {
      af::shared<double> total_intensity(data.size());
      af::const_ref<Shoebox<>> shoebox = data["shoebox"];
      af::shared<bool> success = data["summation_success"];
      af::shared<int> nbg = data["num_pixels.background"];
      af::shared<int> bg_used = data["num_pixels.background_used"];
      af::shared<int> foreground = data["num_pixels.foreground"];
      af::shared<int> valid = data["num_pixels.valid"];
      af::shared<double> variance = data["intensity.sum.variance"];
      af::shared<double> background = data["background.mean"];
      af::shared<double> background_total =
        data["background.sum.value"];  // this is the sum of the background under the
                                       // foreground region
      af::shared<double> background_variance = data["background.sum.variance"];
      for (int i = 0; i < data.size(); i++) {
        background[i] = shoebox[i].mean_background;
        // background_variance[i] = shoebox[i].mean_background;
        double bg_total = shoebox[i].mean_background * shoebox[i].n_valid_fg;
        total_intensity[i] = shoebox[i].total_intensity - bg_total;
        if (shoebox[i].n_invalid_fg > 0) {
          success[i] = false;
        }
        background_total[i] = bg_total;
        double m_n = shoebox[i].n_valid_bg > 0
                       ? (double)shoebox[i].n_valid_fg / (double)shoebox[i].n_valid_bg
                       : 0.0;
        variance[i] = std::abs(total_intensity[i]) + std::abs(bg_total) * (1.0 + m_n);
        background_variance[i] = std::abs(bg_total) * (1.0 + m_n);
        nbg[i] = shoebox[i].n_valid_bg;
        bg_used[i] = shoebox[i].n_valid_bg;
        foreground[i] = shoebox[i].n_valid_fg;
        valid[i] = shoebox[i].n_valid_bg + shoebox[i].n_valid_fg;
      }
      return total_intensity;
    }

    /** @returns The first frame.  */
    int frame0() const {
      return frame0_;
    }

    /** @returns The last frame */
    int frame1() const {
      return frame1_;
    }

    /** @returns The current frame. */
    int frame() const {
      return frame_;
    }

    /** @returns The number of frames  */
    std::size_t nframes() const {
      return nframes_;
    }

    /** @returns The number of panels */
    std::size_t npanels() const {
      return npanels_;
    }

    /**
     * @returns Is the extraction finished.
     */
    bool finished() const {
      return frame_ == frame1_;
    }

    /**
     * @returns The extract time
     */
    double extract_time() const {
      return extract_time_;
    }

    /**
     * @returns The process time
     */
    double process_time() const {
      return process_time_;
    }

  private:
    /**
     * Get an index array specifying which reflections are recorded on a given
     * frame and panel.
     * @param frame The frame number
     * @param panel The panel number
     * @returns An array of indices
     */
    af::const_ref<std::size_t> indices(int frame, std::size_t panel) const {
      std::size_t j0 = panel + (frame - frame0_) * npanels_;
      DIALS_ASSERT(offset_.size() > 0);
      DIALS_ASSERT(j0 < offset_.size() - 1);
      std::size_t i0 = offset_[j0];
      std::size_t i1 = offset_[j0 + 1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - off;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t>(&indices_[off], num);
    }

    af::reflection_table data_;
    double extract_time_;
    double process_time_;
    bool flatten_;
    bool save_;
    std::size_t npanels_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
    double phi0_;
    double dphi_;
    vec3<double> s0_;
    vec3<double> m2_;
    dxtbx::model::Detector detector_;
    double delta_b_r2;
    double delta_m_r2;
    double index0_;
    double index1_;
    std::map<std::size_t, std::vector<profile_model::gaussian_rs::CoordinateSystem>>
      coordinate_systems_;
    std::map<std::size_t, std::vector<double>> attenuation_lengths_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_PROCESSOR_H
