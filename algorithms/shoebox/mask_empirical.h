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
#include <annlib_adaptbx/ann_adaptor.h>
#include <dials/error.h>
#include <boost/math/special_functions/round.hpp>

namespace dials { namespace algorithms { namespace shoebox {

  using annlib_adaptbx::AnnAdaptor;
  using dials::model::Shoebox;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;

  /**
   * A class to mask foreground/background using an empirical approach pixels
   */
  class MaskEmpirical {
  public:
    /**
     * Initialise the stuff needed to create the mask.
     * @param beam The beam model
     * @param detector The detector model
     */
    MaskEmpirical(const af::reflection_table &reference) : reference_(reference) {}

    /**
     * Set all the foreground/background pixels in the shoebox mask by looking at the
     * nearest bright spots
     * @param table Reflection table with a shoebox array and a bbox array for masking
     */
    void mask(af::reflection_table &table) {
      // Convert the reference xyz mm values to a single array for AnnAdaptor, instead
      // of a vec3 double array
      af::shared<vec3<double> > reference_xyz_mm;
      reference_xyz_mm = reference_.get<vec3<double> >("xyzcal.mm");

      af::shared<double> reference_xyz_mm_values;
      for (std::size_t iter = 0; iter < reference_xyz_mm.size(); ++iter) {
        vec3<double> item = reference_xyz_mm[iter];
        reference_xyz_mm_values.push_back(item[0]);
        reference_xyz_mm_values.push_back(item[1]);
        reference_xyz_mm_values.push_back(item[2]);
      }

      // Convert the query xyz mm values to a single array for AnnAdaptor, instead of a
      // vec3 double array
      af::shared<vec3<double> > query_xyz_mm;
      query_xyz_mm = table.get<vec3<double> >("xyzcal.mm");

      af::shared<double> query_xyz_mm_values;
      for (std::size_t iter = 0; iter < query_xyz_mm.size(); ++iter) {
        vec3<double> item = query_xyz_mm[iter];
        query_xyz_mm_values.push_back(item[0]);
        query_xyz_mm_values.push_back(item[1]);
        query_xyz_mm_values.push_back(item[2]);
      }

      // Calculate the nearest neighbors for the query reflections in the reference set
      int nn_window = 10;
      AnnAdaptor A = AnnAdaptor(reference_xyz_mm_values, 3, nn_window);
      A.query(query_xyz_mm_values);

      af::shared<Shoebox<> > reference_shoeboxes =
        reference_.get<Shoebox<> >("shoebox");
      af::shared<Shoebox<> > table_shoeboxes = table.get<Shoebox<> >("shoebox");

      // Generate the mask for each query reflection
      for (std::size_t query_iter = 0; query_iter < table.size(); ++query_iter) {
        Shoebox<> table_shoebox = table_shoeboxes[query_iter];
        af::ref<int, af::c_grid<3> > table_mask = table_shoebox.mask.ref();
        int6 table_bbox = table_shoebox.bbox;
        int tblx1 = table_bbox[0];
        int tbly1 = table_bbox[2];
        int tblz1 = table_bbox[4];
        int tblx2 = table_bbox[1];
        int tbly2 = table_bbox[3];
        int tblz2 = table_bbox[5];

        // Calculate the midpoint of the query reflection.  Center the union of the
        // reference masks around this point
        int tmid_x = boost::math::iround((tblx2 - tblx1) / 2);
        int tmid_y = boost::math::iround((tbly2 - tbly1) / 2);
        int tmid_z = boost::math::iround((tblz2 - tblz1) / 2);

        // Iterate over the nearest bright neigbors of this query reflection
        for (std::size_t nn_iter = query_iter * nn_window;
             nn_iter < (query_iter * nn_window) + nn_window;
             ++nn_iter) {
          Shoebox<> reference_shoebox = reference_shoeboxes[A.nn[nn_iter]];
          af::ref<int, af::c_grid<3> > reference_mask = reference_shoebox.mask.ref();
          int6 bbox_nn = reference_shoebox.bbox;
          int nnx1 = bbox_nn[0];
          int nny1 = bbox_nn[2];
          int nnz1 = bbox_nn[4];
          int nnx2 = bbox_nn[1];
          int nny2 = bbox_nn[3];
          int nnz2 = bbox_nn[5];

          // Calculate the midpoint of the reference reflection.  Center the union of
          // this reflection's mask on the query reflection around this point
          int rmid_x = boost::math::iround((nnx2 - nnx1) / 2);
          int rmid_y = boost::math::iround((nny2 - nny1) / 2);
          int rmid_z = boost::math::iround((nnz2 - nnz1) / 2);

          // Union this reflection's mask with the query reflection
          for (std::size_t z = 0; z < nnz2 - nnz1; z++) {
            for (std::size_t y = 0; y < nny2 - nny1; y++) {
              for (std::size_t x = 0; x < nnx2 - nnx1; x++) {
                table_mask(
                  tmid_z - rmid_z + z, tmid_y - rmid_y + y, tmid_x - rmid_x + x) |=
                  reference_mask(z, y, x);
              }
            }
          }
        }

        // Finally, need to explicitly set the background flags for the rest of the
        // pixels as they are not set by spotfinder in the reference set
        for (std::size_t z = 0; z < tblz2 - tblz1; z++)
          for (std::size_t y = 0; y < tbly2 - tbly1; y++)
            for (std::size_t x = 0; x < tblx2 - tblx1; x++)
              if ((table_mask(z, y, x) & Foreground) != Foreground)
                table_mask(z, y, x) |= Background;
      }
    }

  private:
    af::reflection_table reference_;
  };

}}}  // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_EMPIRICAL_H */
