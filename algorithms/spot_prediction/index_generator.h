/*
 * index_generator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_INDEX_GENERATOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_INDEX_GENERATOR_H

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/loops.h>
#include <dials/array_family/import_scitbx_af.h>

namespace dials { namespace algorithms {

  /**
   * Generate all the possible reflection indices from a unit cell object and
   * calculate the symmetry reduced set of reflections.
   */
  class IndexGenerator {
  public:
    /** Default constructor */
    IndexGenerator() {}

    /**
     * Initialise the generator.
     * @param unit_cell The unit cell structure
     * @param space_group_type The space group type structure
     * @param d_min The resolution
     */
    IndexGenerator(cctbx::uctbx::unit_cell const& unit_cell,
                   cctbx::sgtbx::space_group_type const& space_group_type,
                   double d_min)
        : unit_cell_(unit_cell),
          space_group_type_(space_group_type),
          d_min_(d_min),
          loop_(initialise_loop(unit_cell.max_miller_indices(d_min))) {}

    /**
     * Get the next miller index from the range of miller valid miller indices.
     * Loop over those miller indices specified by the loop structure, if the
     * current index is within the resolution and not systematically absent,
     * then return the miller index, otherwise iterate until a miller index
     * matching those conditions is found. If none are found, then return
     * the index (0, 0, 0).
     * @returns The next miller index
     */
    cctbx::miller::index<> next() {
      for (; loop_.over() == 0;) {
        cctbx::miller::index<> h = loop_();
        loop_.incr();
        if (unit_cell_.d(h) >= d_min_) {
          if (!space_group_type_.group().is_sys_absent(h)) {
            return h;
          }
        }
      }
      return cctbx::miller::index<>(0, 0, 0);
    }

    /**
     * Create an array of miller indices by calling next until (000) is reached.
     * @returns The array of valid miller indices
     */
    af::shared<cctbx::miller::index<> > to_array() {
      af::shared<cctbx::miller::index<> > result;
      for (;;) {
        cctbx::miller::index<> h = this->next();
        if (h.is_zero()) {
          break;
        }
        result.push_back(h);
      }
      return result;
    }

  private:
    /**
     * Initialise the loop to iterate from -max_hkl -> max_hkl
     * @param reference_h_max The maximum miller indices
     * @returns The loop structure
     */
    af::nested_loop<cctbx::miller::index<> > initialise_loop(
      cctbx::miller::index<> const& reference_h_max) const {
      return af::nested_loop<cctbx::miller::index<> >(-reference_h_max,
                                                      reference_h_max + 1);
    }

  private:
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    double d_min_;
    af::nested_loop<cctbx::miller::index<> > loop_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_INDEX_GENERATOR_H
