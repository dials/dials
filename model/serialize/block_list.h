/*
 * block_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_MODEL_SERIALIZE_BLOCK_LIST_H
#define DIALS_MODEL_SERIALIZE_BLOCK_LIST_H

#include <algorithm>
#include <deque>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/basic_shoebox.h>

namespace dials { namespace model {

  /**
   * Class to help identifiy the blocks and shoeboxes in a given range.
   */
  class BlockListIndex {
  public:

    /**
     * Setup a list of block ranges.
     */
    BlockListIndex(const af::const_ref<int2> &blocks)
      : blocks_(blocks.begin(), blocks.end()) {}

    /**
     * Given the shoeboxes, select those to use
     */
    af::shared<std::size_t> indices_in_range(
        int2 zrange, const af::const_ref<double> &z) const {
      af::shared<std::size_t> result;
      DIALS_ASSERT(zrange[0] < zrange[1]);
      for (std::size_t i = 0; i < z.size(); ++i) {
        if (z[i] >= zrange[0] && z[i] < zrange[1]) {
          result.push_back(i);
        }
      }
      return result;
    }

    /**
     * Get the block indices in the given range.
     */
    af::shared<std::size_t> blocks_in_range(int2 z) const {
      af::shared<std::size_t> result;
      for (std::size_t i = 0; i < blocks_.size(); ++i) {
        if ((z[0] >= blocks_[i][0] && z[0] < blocks_[i][1]) ||
            (z[1] <= blocks_[i][1] && z[1] > blocks_[i][0])) {
          result.push_back(i);
        }
      }
      return result;
    }

  private:

    af::shared<int2> blocks_;
  };

  /**
   * A class to help extract shoeboxes in an efficient manner.
   *
   * Start by sorting all the bounding boxes by lower frame number. Then as each
   * frame is added, allocate those shoeboxes which start on the given frame and
   * sort the allocated (or active) shoeboxes by end frame. Now distribute the
   * image pixels to all the active shoeboxes. Finally return a list of
   * shoeboxes that have been completed with the current frame and erase them
   * from the internal list of shoeboxes.
   */
  class BlockList {
  public:

    typedef std::pair<std::size_t, BasicShoebox> record_type;

    /**
     * The block return type
     */
    struct Block {
      int2 zrange;
      af::shared<std::size_t> index;
      af::shared<BasicShoebox> shoebox;
    };

    /**
     * Sort the shoeboxes by starting z
     */
    struct sort_by_start {
      bool operator()(const record_type &a, const record_type &b) const {
        return a.second.bbox[4] < b.second.bbox[4];
      }
    };

    /**
     * Sort the shoeboxes by ending z
     */
    struct sort_by_end {
      bool operator()(const record_type &a, const record_type &b) const {
        return a.second.bbox[5] < b.second.bbox[5];
      }
    };

    /**
     * Sort the index
     */
    struct sort_by_index {
      bool operator()(const record_type &a, const record_type &b) const {
        return a.first < b.first;
      }
    };

    /**
     * Initialise with a list of panels and bounding boxes
     * @param panel The list of panels
     * @param bbox The list of bounding boxes
     * @param z The starting frame number
     */
    BlockList(
        const af::const_ref<std::size_t> &panel,
        const af::const_ref<int6> &bbox,
        int2 zrange)
      : z0_(zrange[0]),
        z_(zrange[0]) {
      DIALS_ASSERT(zrange[1] > zrange[0]);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        DIALS_ASSERT(bbox[i][4] >= zrange[0] && bbox[i][5] <= zrange[1]);
        sbox_.push_back(record_type(i, BasicShoebox(panel[i], bbox[i])));
      }
      std::sort(sbox_.begin(), sbox_.end(), sort_by_start());
      active_ = 0;
    }

    /**
     * Distribute the next frame data and get a list of shoeboxes that have been
     * completely filled and can be written to disk.
     * @param image A list of images
     * @returns The shoeboxes that are done.
     */
    Block next(const af::const_ref<
        af::const_ref< int, af::c_grid<2> > > &image) {
      allocate();
      distribute(image);
      z_++;
      return pop_block();
    }

    /**
     * @returns Are all the shoeboxes finsihed.
     */
    bool empty() const {
      return sbox_.empty();
    }

    /**
     * @returns The current frame number
     */
    int z() const {
      return z_;
    }

  private:

    void allocate() {
      for (; active_ < sbox_.size(); ++active_) {
        BasicShoebox &sbox = sbox_[active_].second;
        if (sbox.bbox[4] > z_) {
          break;
        }
        sbox.allocate();
      }
      std::sort(sbox_.begin(), sbox_.begin() + active_, sort_by_end());
    }

    void distribute(const af::const_ref<
        af::const_ref< int, af::c_grid<2> > > &image) {
      for (std::size_t i = 0; i < active_; ++i) {
        BasicShoebox sbox = sbox_[i].second;
        DIALS_ASSERT(sbox.panel < image.size());
        distribute(sbox, image[sbox.panel]);
      }
    }

    void distribute(BasicShoebox &shoebox,
        const af::const_ref< int, af::c_grid<2> > &image) {

      // Loop through all the indices for this frame
      int2 image_size = image.accessor();

      // Get a reference to a reflection
      int6 bbox = shoebox.bbox;
      int i0 = bbox[0], i1 = bbox[1];
      int j0 = bbox[2], j1 = bbox[3];
      int k0 = bbox[4], k1 = bbox[5];
      int k = z_ - k0;
      DIALS_ASSERT(k0 <= z_ && z_ < k1);
      DIALS_ASSERT(shoebox.is_consistent());

      // Readjust the area to loop over to ensure we're within image bounds
      int jj0 = j0 >= 0 ? j0 : 0;
      int ii0 = i0 >= 0 ? i0 : 0;
      int jj1 = j1 <= image_size[0] ? j1 : image_size[0];
      int ii1 = i1 <= image_size[1] ? i1 : image_size[1];

      // Copy the image pixels
      af::ref< int, af::c_grid<3> > profile = shoebox.data.ref();
      for (int jj = jj0; jj < jj1; ++jj) {
        for (int ii = ii0; ii < ii1; ++ii) {
          int j = jj - j0;
          int i = ii - i0;
          profile(k, j, i) = image(jj, ii);
        }
      }
    }

    Block pop_block() {
      Block result;
      result.zrange[0] = z0_;
      result.zrange[1] = z_;
      std::size_t finished = 0;
      for (; finished < active_; ++finished) {
        std::size_t index = sbox_[finished].first;
        BasicShoebox &sbox = sbox_[finished].second;
        if (sbox.bbox[5] > z_) {
          break;
        } else {
          DIALS_ASSERT(sbox.bbox[5] == z_);
        }
        if (sbox.bbox[4] < result.zrange[0]) {
          result.zrange[0] = sbox.bbox[4];
        }
        result.index.push_back(index);
        result.shoebox.push_back(sbox);
      }
      sbox_.erase(sbox_.begin(), sbox_.begin() + finished);
      active_ -= finished;
      return result;
    }

    int z0_, z_;
    std::deque<record_type> sbox_;
    std::size_t active_;
  };

}} // namespace dials::model


#endif // DIALS_MODEL_SERIALIZE_BLOCK_LIST_H
