/*
 * shoebox_extractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_SHOEBOX_EXTRACTOR_H
#define DIALS_ARRAY_FAMILY_SHOEBOX_EXTRACTOR_H

#include <omptbx/omp_or_stubs.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>

namespace dials { namespace af {

  using model::Image;
  using model::Shoebox;
  using model::Valid;

  /**
   * A class to extract shoebox pixels from images
   */
  class ShoeboxExtractor {
  public:
    /**
     * Initialise the index array. Determine which reflections are recorded on
     * each frame and panel ahead of time to enable quick lookup of the
     * reflections to be written to when processing each image.
     */
    ShoeboxExtractor(af::reflection_table data,
                     std::size_t npanels,
                     int frame0,
                     int frame1)
        : npanels_(npanels),
          frame0_(frame0),
          frame1_(frame1),
          frame_(frame0),
          nframes_(frame1 - frame0) {
      DIALS_ASSERT(frame0_ < frame1_);
      DIALS_ASSERT(npanels_ > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("panel"));
      DIALS_ASSERT(data.contains("bbox"));
      DIALS_ASSERT(data.contains("shoebox"));
      shoebox_ = data["shoebox"];
      af::const_ref<std::size_t> panel = data["panel"];
      af::const_ref<int6> bbox = data["bbox"];
      std::size_t size = nframes_ * npanels_;
      std::vector<std::size_t> num(size, 0);
      std::vector<std::size_t> count(size, 0);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(bbox[i][4] >= frame0_);
        DIALS_ASSERT(bbox[i][5] <= frame1_);
        DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
        DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
          std::size_t j = panel[i] + (z - frame0_) * npanels_;
          DIALS_ASSERT(j < num.size());
          num[j]++;
        }
      }
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));
      indices_.resize(offset_.back());
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
          std::size_t j = panel[i] + (z - frame0_) * npanels_;
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
    void next(const Image<T>& image) {
      typedef Shoebox<>::float_type float_type;
      typedef af::ref<float_type, af::c_grid<3> > sbox_data_type;
      typedef af::ref<int, af::c_grid<3> > sbox_mask_type;
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == npanels_);
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref<T, af::c_grid<2> > data = image.data(p);
        af::const_ref<bool, af::c_grid<2> > mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox_.size());
          Shoebox<>& sbox = shoebox_[ind[i]];
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
          std::size_t xs = x1 - x0;
          std::size_t ys = y1 - y0;
          std::size_t z = frame_ - z0;
          std::size_t yi = data.accessor()[0];
          std::size_t xi = data.accessor()[1];
          int xb = x0 >= 0 ? 0 : std::abs(x0);
          int yb = y0 >= 0 ? 0 : std::abs(y0);
          int xe = x1 <= xi ? xs : xs - (x1 - (int)xi);
          int ye = y1 <= yi ? ys : ys - (y1 - (int)yi);
          DIALS_ASSERT(ye > yb && yb >= 0 && ye <= ys);
          DIALS_ASSERT(xe > xb && xb >= 0 && xe <= xs);
          DIALS_ASSERT(yb + y0 >= 0 && ye + y0 <= yi);
          DIALS_ASSERT(xb + x0 >= 0 && xe + x0 <= xi);
          DIALS_ASSERT(sbox.is_consistent());
          for (std::size_t y = yb; y < ye; ++y) {
            for (std::size_t x = xb; x < xe; ++x) {
              sdata(z, y, x) = data(y + y0, x + x0);
              smask(z, y, x) = mask(y + y0, x + x0) ? Valid : 0;
            }
          }
        }
      }
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

    std::size_t npanels_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    af::shared<Shoebox<> > shoebox_;
    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
  };

}}  // namespace dials::af

#endif  // DIALS_ARRAY_FAMILY_SHOEBOX_EXTRACTOR_H
