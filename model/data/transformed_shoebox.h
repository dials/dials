/*
 * transformed_shoebox.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_TRANSFORMED_SHOEBOX_H
#define DIALS_MODEL_DATA_TRANSFORMED_SHOEBOX_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace model {

  /**
   * A class to hold transformed shoebox information
   */
  struct TransformedShoebox {

    af::versa< double, af::c_grid<3> > data;       ///< The shoebox data
    af::versa< double, af::c_grid<3> > background; ///< The shoebox background

    /**
     * Initialise the shoebox
     */
    TransformedShoebox()
      : data(af::c_grid<3>(0, 0, 0)),
        background(af::c_grid<3>(0, 0, 0)) {}

    /** @return True/False whether the array and bbox sizes are consistent */
    bool is_consistent() const {
      return data.accessor().all_eq(background.accessor());
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are the same
     */
    bool operator==(const TransformedShoebox &rhs) const {
      return data.all_eq(rhs.data) &&
             background.all_eq(rhs.background);
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are not the same
     */
    bool operator!=(const TransformedShoebox &rhs) const {
      return !(*this == rhs);
    }
  };

}};

#endif /* DIALS_MODEL_DATA_TRANSFORMED_SHOEBOX_H */
