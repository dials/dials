/*
 * populator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H
#define DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H

#include <boost/unordered_map.hpp>
#include <scitbx/array_family/tiny_types.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using boost::unordered_map;
  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::int6;
  using dials::model::Shoebox;
  using dials::model::Valid;

  class Extractor {
  public:

    Extractor(const af::const_ref<std::size_t> &panels,
              const af::const_ref<int6> &bboxes,
              const af::shared<bool> &mask,
              const af::shared<double> &gain,
              const af::shared<double> &dark,
              const af::shared<int2> &shape)
      : shoeboxes_(panels.size()),
        gain_(gain),
        dark_(dark),
        shape_(shape),
        offset_(shape.size()),
        npanels_(shape.size()) {
      DIALS_ASSERT(panels.size() == bboxes.size());
      for (std::size_t i = 0; i < shoeboxes_.size(); ++i) {
        shoeboxes_[i] = Shoebox(panels[i], bboxes[i]);
        shoeboxes_[i].allocate();
      }
      DIALS_ASSERT(shape_.size() > 0);
      offset_[0] = 0;
      for (std::size_t i = 0; i < shape_.size(); ++i) {
        DIALS_ASSERT(shape_[i].all_gt(0));
        offset_[i+1] = offset_[i] + af::product(shape_[i]);
      }
      DIALS_ASSERT(mask.size() == offset_.back() + af::product(shape_.back()));
      DIALS_ASSERT(mask.size() == gain_.size());
      DIALS_ASSERT(mask.size() == dark_.size());
      initialise(mask.const_ref());
    }

    Extractor(const af::const_ref<int6> &bboxes,
              af::versa<bool, af::c_grid<2> > mask,
              af::versa<double, af::c_grid<2> > gain,
              af::versa<double, af::c_grid<2> > dark)
      : shoeboxes_(bboxes.size()),
        gain_(gain.as_1d()),
        dark_(dark.as_1d()),
        shape_(1, mask.accessor()),
        offset_(1, 0),
        npanels_(1) {
      for (std::size_t i = 0; i < shoeboxes_.size(); ++i) {
        shoeboxes_[i] = Shoebox(bboxes[i]);
        shoeboxes_[i].allocate();
      }
      DIALS_ASSERT(mask.accessor().all_gt(0));
      DIALS_ASSERT(mask.accessor().all_eq(gain.size()));
      DIALS_ASSERT(mask.accessor().all_eq(dark.size()));
      initialise(mask.const_ref().as_1d());
    }

    void add_image(std::size_t panel, std::size_t frame,
      const af::const_ref< int, af::c_grid<2> > &image) {

      // Get the image size
      int2 image_size = image.accessor();
      DIALS_ASSERT(panel >= 0 && panel < npanels_);
      DIALS_ASSERT(image_size.all_eq(shape_[panel]));

      // Get the indices for this frame
      af::shared<int> ind = indices(panel, frame);

      // Loop through all the indices for this frame
      for (std::size_t i = 0; i < ind.size(); ++i) {

        // Get a reference to a reflection
        Shoebox &shoebox = shoeboxes_[ind[i]];
        int6 bbox = shoebox.bbox;
        int i0 = bbox[0], i1 = bbox[1];
        int j0 = bbox[2], j1 = bbox[3];
        int k0 = bbox[4], k1 = bbox[5];
        int k = frame - k0;

        // Readjust the area to loop over to ensure we're within image bounds
        int jj0 = j0 >= 0 ? j0 : 0;
        int ii0 = i0 >= 0 ? i0 : 0;
        int jj1 = j1 <= image_size[0] ? j1 : image_size[0];
        int ii1 = i1 <= image_size[1] ? i1 : image_size[1];

        // Get the reflection profile
        af::ref< double, af::c_grid<3> > profile = shoebox.data.ref();
        DIALS_ASSERT(profile.accessor()[0] == (k1 - k0));
        DIALS_ASSERT(profile.accessor()[1] == (j1 - j0));
        DIALS_ASSERT(profile.accessor()[2] == (i1 - i0));

        // Copy the image pixels
        for (int jj = jj0; jj < jj1; ++jj) {
          for (int ii = ii0; ii < ii1; ++ii) {
            int j = jj - j0;
            int i = ii - i0;
            int kk = offset_[panel] + ii + jj * image_size[1];
            profile(k, j, i) = gain_[kk] * (image(jj, ii) - dark_[kk]);
          }
        }
      }
    }

    af::shared<Shoebox> shoeboxes() {
      return shoeboxes_;
    }

    af::shared<int> indices(std::size_t panel, std::size_t frame) {
      DIALS_ASSERT(panel >= 0 && panel < npanels_);
      return indices_[panel + frame * npanels_];
    }

  private:


    void initialise(const af::const_ref<bool> &detector_mask) {
      initialise_indices();
      initialise_mask(detector_mask);
    }

    void initialise_indices() {
      for (std::size_t i = 0; i < shoeboxes_.size(); ++i) {
        Shoebox& shoebox = shoeboxes_[i];
        DIALS_ASSERT(shoebox.panel >= 0 && shoebox.panel < npanels_);
        for (int f = shoebox.bbox[4]; f < shoebox.bbox[5]; ++f) {
          indices_[shoebox.panel + f * npanels_].push_back(i);
        }
      }
    }

    void initialise_mask(const af::const_ref<bool> &detector_mask) {

      // Loop through the shoeboxes
      for (std::size_t s = 0; s < shoeboxes_.size(); ++s) {

        // Loop through all the points in the shoebox mask. If the point is
        // within the detector mask area, get the detector mask value and
        // set the shoebox mask value to Valid, otherwise set the shoebox
        // mask value to 0.
        af::ref< int, af::c_grid<3> > mask = shoeboxes_[s].mask.ref();
        int6 bbox = shoeboxes_[s].bbox;
        std::size_t panel = shoeboxes_[s].panel;
        int2 isize = shape_[panel];
        for (std::size_t j = 0; j < mask.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < mask.accessor()[2]; ++i) {
            int di = bbox[0] + i, dj = bbox[2] + j;
            bool mask_value = false;
            if (di >= 0 && di < isize[1] && dj >= 0 && dj < isize[0]) {
              mask_value = detector_mask[offset_[panel] + di + dj * isize[1]];
            }
            for (std::size_t k = 0; k < mask.accessor()[0]; ++k) {
              mask(k, j, i) = (mask_value ? Valid : 0);
            }
          }
        }
      }
    }

    af::shared<Shoebox> shoeboxes_;
    af::shared<double> gain_;
    af::shared<double> dark_;
    af::shared<int2> shape_;
    af::shared<int> offset_;
    std::size_t npanels_;
    unordered_map< int, af::shared<int> > indices_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_POPULATOR_H */
