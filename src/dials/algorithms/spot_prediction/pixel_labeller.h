/*
 * pixel_labeller.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_LABELLER_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_LABELLER_H

#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * A class to do pixel labelling as in XDS.
   */
  class PixelLabeller {
  public:
    typedef af::versa<vec3<double>, af::c_grid<2> > pstar_array_type;
    typedef af::shared<pstar_array_type> array_type;

    /**
     * Preprocess the beam and detector
     * @param beam The beam model
     * @param detector The detector model
     */
    PixelLabeller(BeamBase &beam, Detector detector) {
      p_star_.resize(detector.size());
      vec3<double> s0 = beam.get_s0();
      for (std::size_t p = 0; p < detector.size(); ++p) {
        const Panel &panel = detector[p];
        vec2<std::size_t> image_size = panel.get_image_size();
        p_star_[p].resize(af::c_grid<2>(image_size[1], image_size[0]));
        pstar_array_type ps = p_star_[p];
        for (std::size_t j = 0; j < image_size[1]; ++j) {
          for (std::size_t i = 0; i < image_size[0]; ++i) {
            vec3<double> s1 = panel.get_pixel_lab_coord(vec2<double>(j + 0.5, i + 0.5));
            ps(j, i) = s1 - s0;
          }
        }
      }
    }

    /**
     * @returns The number of panels.
     */
    std::size_t size() const {
      return p_star_.size();
    }

    /**
     * @returns The size of the requested panel.
     */
    af::c_grid<2> panel_size(std::size_t index) const {
      DIALS_ASSERT(index < size());
      return p_star_[index].accessor();
    }

    /**
     * Label the pixels using the given crystal A matrix.
     * @param index The index labels
     * @param A The setting matrix
     * @param panel_number The panel
     */
    void label(af::ref<cctbx::miller::index<> > index,
               mat3<double> A,
               std::size_t panel_number) const {
      DIALS_ASSERT(panel_number < size());
      af::c_grid<2> size = panel_size(panel_number);
      DIALS_ASSERT(index.size() == size[0] * size[1]);
      mat3<double> A1 = A.inverse();
      pstar_array_type ps = p_star_[panel_number];
      for (std::size_t j = 0; j < size[0]; ++j) {
        for (std::size_t i = 0; i < size[1]; ++i) {
          vec3<double> hf = A1 * ps(j, i);
          cctbx::miller::index<> h((int)std::floor(hf[0] + 0.5),
                                   (int)std::floor(hf[1] + 0.5),
                                   (int)std::floor(hf[2] + 0.5));
          index[i + j * size[1]] = h;
        }
      }
    }

  private:
    array_type p_star_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_PIXEL_LABELLER_H
