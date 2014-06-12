/*
 * mask_empirical.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_MASK_EMPIRICAL_H
#define DIALS_ALGORITHMS_SHOEBOX_MASK_EMPIRICAL_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/array_family/reflection_table.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using dials::model::Shoebox;

  /**
   * A class to mask foreground/background using an empirical approach pixels
   */
  class MaskEmpirical{
  public:

    /**
     * Initialise the stuff needed to create the mask.
     * @param beam The beam model
     * @param detector The detector model
     */
    MaskEmpirical(const af::reflection_table &reference)
      : reference_(reference){
    }

    /**
     * Set all the foreground/background pixels in the shoebox mask.
     * @param table Reflection table with a shoebox array and a bbox array for masking
     */
    void mask(af::reflection_table &table) {
      printf("Not implemented\n");
    }


  private:
    af::reflection_table reference_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_EMPIRICAL_H */
