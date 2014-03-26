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
  private:

    /**
     * Sort the shoeboxes by starting z
     */
    struct sort_by_start {
      bool operator()(const BasicShoebox &a, const BasicShoebox &b) const {
        return a.bbox[4] < b.bbox[4];
      }
    };

    /**
     * Sort the shoeboxes by ending z
     */
    struct sort_by_end {
      bool operator()(const BasicShoebox &a, const BasicShoebox &b) const {
        return a.bbox[5] < b.bbox[5];
      }
    };

  public:

    /**
     * Initialise with a list of panels and bounding boxes
     * @param panel The list of panels
     * @param bbox The list of bounding boxes
     * @param z The starting frame number
     */
    BlockList(
        const af::const_ref<std::size_t> &panel,
        const af::const_ref<int6> &bbox, int z)
      : z_(z) {
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        sbox_.push_back(BasicShoebox(i, panel[i], bbox[i]));
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
    af::shared<BasicShoebox> next(
        const af::const_ref< af::const_ref< int, af::c_grid<2> > > &image) {
      allocate();
      distribute(image);
      af::shared<BasicShoebox> result = pop_block();
      z_++;
      return result;
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
        BasicShoebox &sbox = sbox_[active_];
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
        std::size_t p = sbox_[i].panel;
        DIALS_ASSERT(p < image.size());
        distribute(sbox_[i], image[p]);
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

    af::shared<BasicShoebox> pop_block() {
      af::shared<BasicShoebox> result;
      for (std::size_t i = 0; i < active_; ++i) {
        BasicShoebox &sbox = sbox_[i];
        if (sbox.bbox[5] > z_) {
          break;
        } else {
          DIALS_ASSERT(sbox.bbox[5] == z_);
        }
        result.push_back(sbox);
        sbox_.pop_front();
        --active_;
      }
      return result;
    }

    int z_;
    std::deque<BasicShoebox> sbox_;
    std::size_t active_;
  };

}} // namespace dials::model


#endif // DIALS_MODEL_SERIALIZE_BLOCK_LIST_H
