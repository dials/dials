/*
 * model.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_GLM_MODEL_H
#define DIALS_ALGORITHMS_BACKGROUND_GLM_MODEL_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * Base class for background models
   */
  class BackgroundModel {
  public:
    virtual ~BackgroundModel(){};

    virtual af::versa<double, af::c_grid<3> > extract(std::size_t panel,
                                                      int6 bbox) const = 0;
  };

  /**
   * A simple static background model
   */
  class StaticBackgroundModel : public BackgroundModel {
  public:
    StaticBackgroundModel() {}

    /**
     * Extract a shoebox
     * @param bbox The bounding box
     * @returns The model data
     */
    virtual af::versa<double, af::c_grid<3> > extract(std::size_t panel,
                                                      int6 bbox) const {
      DIALS_ASSERT(panel < data_.size());
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
      af::c_grid<3> grid(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);
      af::versa<double, af::c_grid<3> > result(grid, 0);
      af::const_ref<double, af::c_grid<2> > data = data_[panel].const_ref();
      for (std::size_t j = 0; j < result.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < result.accessor()[2]; ++i) {
          int ii = bbox[0] + i;
          int jj = bbox[2] + j;
          if (ii >= 0 && jj >= 0 && ii < data.accessor()[1]
              && jj < data.accessor()[0]) {
            double value = data(jj, ii);
            for (std::size_t k = 0; k < result.accessor()[0]; ++k) {
              result(k, j, i) = value;
            }
          }
        }
      }
      return result;
    }

    /**
     * Add the background model
     * @param data The model data
     */
    void add(const af::const_ref<double, af::c_grid<2> > &data) {
      af::versa<double, af::c_grid<2> > temp(data.accessor());
      std::copy(data.begin(), data.end(), temp.begin());
      data_.push_back(temp);
    }

    /**
     * The number of panels
     */
    std::size_t size() const {
      return data_.size();
    }

    /**
     * Get the data array
     * @returns The data array
     */
    af::versa<double, af::c_grid<2> > data(std::size_t panel) const {
      DIALS_ASSERT(panel < size());
      return data_[panel];
    }

  protected:
    af::shared<af::versa<double, af::c_grid<2> > > data_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_BACKGROUND_GLM_MODEL_H
