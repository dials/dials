/*
 * shoebox_block_importer.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_SERIALIZE_SHOEBOX_BLOCK_IMPORTER_H
#define DIALS_MODEL_SERIALIZE_SHOEBOX_BLOCK_IMPORTER_H

#include <vector>
#include <dials/model/serialize/shoebox_file_importer.h>

namespace dials { namespace model { namespace serialize {

  /**
   * A class to import shoeboxes in blocks of frames
   */
  class ShoeboxBlockImporter {
  public:

    typedef ShoeboxFileImporter::gain_map_ref_type gain_map_ref_type;
    typedef ShoeboxFileImporter::dark_map_ref_type dark_map_ref_type;
    typedef ShoeboxFileImporter::mask_map_ref_type mask_map_ref_type;

    /**
     * Initialise the importer
     * @param filename The path to the shoebox file
     */
    ShoeboxBlockImporter(const std::string &filename,
                         const af::const_ref<std::size_t> &blocks)
        : importer_(filename) {
      init_indices(blocks);
    }

    /**
     * Initialise the importer
     * @param filename The path to the shoebox file
     * @param gain The gain map
     * @param dark The dark map
     * @param mask The mask
     */
    ShoeboxBlockImporter(const std::string &filename,
                         const af::const_ref<std::size_t> &blocks,
                         const af::const_ref<gain_map_ref_type> &gain,
                         const af::const_ref<dark_map_ref_type> &dark,
                         const af::const_ref<mask_map_ref_type> &mask)
        : importer_(filename, gain, dark, mask) {
      init_indices(blocks);
    }

    /**
     * @returns The number of shoeboxes
     */
    std::size_t size() const {
      return indices_.size();
    }

    /**
     * Read a shoebox
     * @param index The index of the shoebox
     * @returns The indices and shoebox list
     */
    std::pair< af::shared<std::size_t>, af::shared< Shoebox<> > >
    operator[](std::size_t index) {
      DIALS_ASSERT(index < size());
      af::shared<std::size_t> ind = indices_[index];
      return std::make_pair(
          af::shared<std::size_t>(ind.begin(), ind.end()),
          importer_.select(ind.const_ref()));
    }

    /**
     * Read a blob of data
     */
    std::string blob() {
      return importer_.blob();
    }

  private:

    void init_indices(const af::const_ref<std::size_t> &blocks) {

      // Ensure blocks are in order
      DIALS_ASSERT(blocks.size() >= 2);
      for (std::size_t i = 1; i < blocks.size(); ++i) {
        DIALS_ASSERT(blocks[i] > blocks[i-1]);
      }

      // Ensure centroids are in order
      af::shared<double> z = importer_.z();
      for (std::size_t i = 1; i < z.size(); ++i) {
        DIALS_ASSERT(z[i] >= z[i-1]);
      }

      // Create the list of indices
      indices_.reserve(blocks.size() - 1);
      std::size_t j = 0;
      for (; j < z.size() && z[j] < blocks[0]; ++j);
      for (std::size_t i = 0; i < blocks.size() - 1; ++i) {
        af::shared<std::size_t> ind;
        for (; j < z.size() && z[j] < blocks[i+1]; ++j) {
          ind.push_back(j);
        }
        indices_.push_back(ind);
      }
    }

    ShoeboxFileImporter importer_;
    af::shared< af::shared<std::size_t> > indices_;
  };

}}} // namespace dials::model::serialize

#endif // DIALS_MODEL_SERIALIZE_SHOEBOX_BLOCK_IMPORTER_H
