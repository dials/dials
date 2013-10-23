/*
 * partial_extractor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_PARTIAL_EXTRACTOR_H
#define DIALS_ALGORITHMS_SHOEBOX_PARTIAL_EXTRACTOR_H

#include <vector>
#include <scitbx/array_family/tiny_types.h>
#include <dials/model/data/partial_shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::int6;
  using dials::model::PartialShoebox;

  /**
   * A class to extract shoeboxes from a sweep
   */
  class PartialExtractor {
  public:

    /**
     * Initialise the shoeboxes with the given input panels and bounding boxes
     * Setup the shoebox indices and the mask with the given detector mask.
     * This constructor is for multiple panels.
     * @param panels The list of panel numbers for the shoeboxes
     * @param bboxes The list of bounding boxes for the shoeboxes
     * @param zrange The range of frames to extract
     * @param npanels The number of panels
     */
    PartialExtractor(const af::const_ref<std::size_t> &panels,
                     const af::const_ref<int6> &bboxes,
                     int2 zrange,
                     std::size_t npanels)
      : npanels_(npanels),
        zrange_(zrange) {
      DIALS_ASSERT(zrange[0] < zrange[1]);
      DIALS_ASSERT(bboxes.size() > 0);
      DIALS_ASSERT(panels.size() == bboxes.size());
      indices_.resize(af::c_grid<2>(npanels_, zrange_[1] - zrange_[0]));
      for (std::size_t i = 0; i < bboxes.size(); ++i) {
        DIALS_ASSERT(panels[i] < npanels_);
        int z0 = std::max(zrange[0], bboxes[i][4]);
        int z1 = std::min(zrange[1], bboxes[i][5]);
        if (z0 < z1) {
          PartialShoebox s(panels[i], bboxes[i], int2(z0, z1));
          s.allocate();
          DIALS_ASSERT(s.is_consistent());
          shoeboxes_.push_back(s);
          shoebox_indices_.push_back(i);
        }
      }
      initialise_indices();
    }

    /**
     * Initialise the shoeboxes with the given input bounding boxes
     * Setup the shoebox indices and the mask with the given detector mask.
     * This constructor is for single panels.
     * @param bboxes The list of bounding boxes for the shoeboxes
     * @param zrange The number of frames to extract
     */
    PartialExtractor(const af::const_ref<int6> &bboxes, int2 zrange)
      : npanels_(1),
        zrange_(zrange) {
      DIALS_ASSERT(zrange[0] < zrange[1]);
      DIALS_ASSERT(bboxes.size() > 0);
      indices_.resize(af::c_grid<2>(npanels_, zrange_[1] - zrange_[0]));
      for (std::size_t i = 0; i < bboxes.size(); ++i) {
        int z0 = std::max(zrange[0], bboxes[i][4]);
        int z1 = std::min(zrange[1], bboxes[i][5]);
        if (z0 < z1) {
          PartialShoebox s(bboxes[i], int2(z0, z1));
          s.allocate();
          DIALS_ASSERT(s.is_consistent());
          shoeboxes_.push_back(s);
          shoebox_indices_.push_back(i);
        }
      }
      initialise_indices();
    }

    /**
     * Add the data from a single image.
     * @param panel The panel number
     * @param frame The frame number
     * @param image The image pixels
     */
    void add_image(std::size_t panel, int frame,
      const af::const_ref< int, af::c_grid<2> > &image) {

      // Check the input
      DIALS_ASSERT(panel >= 0 && panel < npanels_);
      DIALS_ASSERT(frame >= zrange_[0] && frame < zrange_[1]);

      // Loop through all the indices for this frame
      int2 image_size = image.accessor();
      const std::vector<int> &ind = indices_(panel, frame - zrange_[0]);
      for (std::size_t i = 0; i < ind.size(); ++i) {

        // Get a reference to a reflection
        PartialShoebox &shoebox = shoeboxes_[ind[i]];
        int6 bbox = shoebox.bbox;
        int2 zrange = shoebox.zrange;
        int i0 = bbox[0], i1 = bbox[1];
        int j0 = bbox[2], j1 = bbox[3];
        int k0 = zrange[0];
        int k = frame - k0;

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
    }

    /** @returns The indices of the shoeboxes into the original arrays */
    af::shared<int> shoebox_indices() const {
      af::shared<int> result(shoebox_indices_.size(), 0);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = shoebox_indices_[i];
      }
      return result;
    }

    /** @returns The list of shoeboxes */
    af::shared<PartialShoebox> shoeboxes() const {
      af::shared<PartialShoebox> result(
        shoeboxes_.size(), PartialShoebox());
      for (std::size_t i = 0; i < shoeboxes_.size(); ++i) {
        result[i] = shoeboxes_[i];
      }
      return result;
    }

    /**
     * @param panel The panel number
     * @param frame The frame number
     * @returns The shoebox indices for the given panel and frame
     */
    af::shared<int> indices(std::size_t panel, int frame) const {
      DIALS_ASSERT(panel >= 0 && panel < npanels_);
      DIALS_ASSERT(frame >= zrange_[0] && frame < zrange_[1]);
      const std::vector<int> &ind = indices_(panel, frame - zrange_[0]);
      af::shared<int> result(ind.size(), af::init_functor_null<int>());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = ind[i];
      }
      return result;
    }

  private:

    /**
     * Initialise the indices of the shoeboxes on each frame
     */
    void initialise_indices() {
      for (std::size_t i = 0; i < shoeboxes_.size(); ++i) {
        PartialShoebox& shoebox = shoeboxes_[i];
        for (int f = shoebox.zrange[0]; f < shoebox.zrange[1]; ++f) {
          indices_(shoebox.panel, f - zrange_[0]).push_back(i);
        }
      }
    }

    std::size_t npanels_;
    int2 zrange_;
    af::shared<PartialShoebox> shoeboxes_;
    af::shared<int> shoebox_indices_;
    af::versa< std::vector<int>, af::c_grid<2> > indices_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_PARTIAL_EXTRACTOR_H */
