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

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <ctime>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <dials/algorithms/integration/compressed_slice_store.h>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dxtbx/array_family/flex_table_suite.h>
#include <dials/array_family/boost_python/reflection_table_suite.h>

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
     *
     * If compress is set, the slices of each 3D shoebox are kept compressed
     * until its last frame, and only then decoded into the shoebox arrays.
     * Flat shoeboxes are a single slice, so they are never compressed.
     */
    ShoeboxProcessor(af::reflection_table data,
                     std::size_t npanels,
                     int frame0,
                     int frame1,
                     bool save,
                     bool compress = false)
        : data_(data),
          extract_time_(0.0),
          process_time_(0.0),
          save_(save),
          compress_(compress),
          batch_bytes_(0),
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
      af::const_ref<Shoebox<> > shoebox = data["shoebox"];
      std::size_t size = nframes_ * npanels_;
      std::vector<std::size_t> num(size, 0);
      std::vector<std::size_t> count(size, 0);
      flatten_ = shoebox[0].flat;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        DIALS_ASSERT(shoebox[i].flat == flatten_);
        DIALS_ASSERT(shoebox[i].is_data_allocated() == false);
        DIALS_ASSERT(shoebox[i].is_background_allocated() == false);
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
      compress_ = compress_ && !flatten_;
      if (compress_) {
        store_ = boost::make_shared<CompressedSliceStore>(shoebox.size());
        finishing_.assign(nframes_, 0);
        for (std::size_t i = 0; i < shoebox.size(); ++i) {
          int z1 = shoebox[i].bbox[5];
          if (z1 > frame0_ && z1 <= frame1_) {
            finishing_[z1 - 1 - frame0_] += decoded_shoebox_bytes(
              shoebox[i].xsize() * shoebox[i].ysize() * shoebox[i].zsize());
          }
        }
        batch_bytes_ = finishing_batch_bytes(finishing_);
      }
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
      typedef af::ref<float_type, af::c_grid<3> > sbox_data_type;
      typedef af::ref<uint8_t, af::c_grid<3> > sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);

      // Get the initial time
      double start_time = timestamp();

      // For each image, extract shoeboxes of reflections recorded.
      // Allocate data where necessary
      af::ref<Shoebox<> > shoebox = data_["shoebox"];
      af::shared<std::size_t> process_indices;
      // With compression, if the shoeboxes finishing on this frame would take
      // too much memory decoded all at once, store their last slices too, and
      // decode and process them in batches afterwards.
      bool defer = compress_ && finishing_[frame_ - frame0_] > batch_bytes_;
      af::shared<std::size_t> unprocessed_indices;
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2> > data = image.data(p);
        af::const_ref<bool, af::c_grid<2> > mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox.size());
          Shoebox<>& sbox = shoebox[ind[i]];
          bool last = frame_ == sbox.bbox[5] - 1;
          bool packed = compress_ && (!last || defer);
          if (frame_ == sbox.bbox[4]) {
            DIALS_ASSERT(sbox.is_data_allocated() == false);
            DIALS_ASSERT(sbox.is_background_allocated() == false);
            if (!packed) {
              sbox.allocate_data();
            }
          } else if (compress_ && last && !packed) {
            // Decode the earlier slices, then write this one in place
            decode(ind[i], frame_ - sbox.bbox[4]);
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
          if (packed) {
            std::size_t n = (std::size_t)xs * ys;
            if (xb == 0 && yb == 0 && xe == xs && ye == ys) {
              // Every pixel is overwritten below
              slice_data_.resize(n);
              slice_mask_.resize(n);
            } else {
              // Pixels off the image stay zero, as in a newly allocated shoebox
              slice_data_.assign(n, 0);
              slice_mask_.assign(n, 0);
            }
          }
          if (yb >= ye || xb >= xe) {
            if (packed) {
              store_->append(
                ind[i], slice_data_.data(), slice_mask_.data(), slice_data_.size());
              if (last) {
                unprocessed_indices.push_back(ind[i]);
              }
            }
            continue;
          }
          DIALS_ASSERT(yb >= 0 && ye <= ys);
          DIALS_ASSERT(xb >= 0 && xe <= xs);
          DIALS_ASSERT(yb + y0 >= 0 && ye + y0 <= yi);
          DIALS_ASSERT(xb + x0 >= 0 && xe + x0 <= xi);
          if (packed) {
            for (std::size_t y = yb; y < ye; ++y) {
              float_type* srow = slice_data_.data() + y * xs;
              uint8_t* mrow = slice_mask_.data() + y * xs;
              for (std::size_t x = xb; x < xe; ++x) {
                srow[x] = data(y + y0, x + x0);
                mrow[x] = mask(y + y0, x + x0) ? Valid : 0;
              }
            }
            store_->append(
              ind[i], slice_data_.data(), slice_mask_.data(), slice_data_.size());
          } else if (flatten_) {
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
      if (defer) {
        // A shoebox whose last slice is off the image is not processed, but
        // is left allocated, as without compression
        for (std::size_t i = 0; i < unprocessed_indices.size(); ++i) {
          const Shoebox<>& sbox = shoebox[unprocessed_indices[i]];
          decode(unprocessed_indices[i], sbox.bbox[5] - sbox.bbox[4]);
        }
        std::size_t i0 = 0;
        while (i0 < process_indices.size()) {
          std::size_t i1 = i0;
          std::size_t bytes = 0;
          while (i1 < process_indices.size()) {
            const Shoebox<>& sbox = shoebox[process_indices[i1]];
            std::size_t next_bytes =
              decoded_shoebox_bytes(sbox.xsize() * sbox.ysize() * sbox.zsize());
            if (i1 > i0 && bytes + next_bytes > batch_bytes_) {
              break;
            }
            bytes += next_bytes;
            ++i1;
          }
          double decode_start = timestamp();
          for (std::size_t i = i0; i < i1; ++i) {
            const Shoebox<>& sbox = shoebox[process_indices[i]];
            decode(process_indices[i], sbox.bbox[5] - sbox.bbox[4]);
          }
          extract_time_ += timestamp() - decode_start;
          process(af::const_ref<std::size_t>(&process_indices[i0], i1 - i0), executor);
          i0 = i1;
        }
      } else if (process_indices.size() > 0) {
        process(process_indices.const_ref(), executor);
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
      typedef af::ref<float_type, af::c_grid<3> > sbox_data_type;
      typedef af::ref<uint8_t, af::c_grid<3> > sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);

      // Get the initial time
      double start_time = timestamp();

      // For each image, extract shoeboxes of reflections recorded.
      // Allocate data where necessary
      af::ref<Shoebox<> > shoebox = data_["shoebox"];
      af::shared<std::size_t> process_indices;
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2> > data = image.data(p);
        af::const_ref<bool, af::c_grid<2> > mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox.size());
          Shoebox<>& sbox = shoebox[ind[i]];
          if (frame_ == sbox.bbox[4]) {
            DIALS_ASSERT(sbox.is_data_allocated() == false);
            DIALS_ASSERT(sbox.is_background_allocated() == false);
            sbox.allocate_data();
            sbox.allocate_background();
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

    /**
     * @returns Whether the shoeboxes are compressed while being read
     */
    bool compressed() const {
      return compress_;
    }

    /**
     * @returns The most memory held at once by compressed slices, in bytes
     */
    std::size_t max_compressed_bytes() const {
      return store_ ? store_->max_bytes_used() : 0;
    }

  private:
    /**
     * Allocate a shoebox's data and mask, and decode its stored slices into
     * them
     * @param index The reflection index
     * @param nslices The number of slices stored
     */
    void decode(std::size_t index, std::size_t nslices) {
      af::ref<Shoebox<> > shoebox = data_["shoebox"];
      Shoebox<>& sbox = shoebox[index];
      sbox.allocate_data();
      std::size_t n = sbox.xsize() * sbox.ysize();
      store_->extract(index, sbox.data.begin(), sbox.mask.begin(), n, nslices);
    }

    /**
     * Run the executor on a batch of complete shoeboxes, then free them
     * @param ind The reflection indices
     * @param executor The executor
     */
    void process(af::const_ref<std::size_t> ind, Executor& executor) {
      using dials::af::boost_python::reflection_table_suite::select_rows_index;
      using dxtbx::af::flex_table_suite::set_selected_rows_index;
      double start_time = timestamp();
      af::ref<Shoebox<> > shoebox = data_["shoebox"];
      // allocate background only for the selected reflections
      for (std::size_t i = 0; i < ind.size(); ++i) {
        shoebox[ind[i]].allocate_background();
      }
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
      return af::const_ref<std::size_t>(indices_.data() + off, num);
    }

    af::reflection_table data_;
    double extract_time_;
    double process_time_;
    bool flatten_;
    bool save_;
    bool compress_;
    std::size_t npanels_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
    boost::shared_ptr<CompressedSliceStore> store_;
    std::vector<std::size_t> finishing_;
    std::size_t batch_bytes_;
    std::vector<Shoebox<>::float_type> slice_data_;
    std::vector<uint8_t> slice_mask_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_PROCESSOR_H
