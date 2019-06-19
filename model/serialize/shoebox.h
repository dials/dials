

#ifndef DIALS_MODEL_SERIALIZE_SHOEBOX_H
#define DIALS_MODEL_SERIALIZE_SHOEBOX_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <string>

namespace dials { namespace model {

  class ShoeboxReaderInternal {
  public:
    ShoeboxReaderInternal(const std::string &filename) {}

    af::shared<Shoebox<double> > read(std::size_t z0, std::size_t z1) const {}

  private:
  };

  /**
   * Interface for reading shoeboxes from the intermediate file.
   */
  class ShoeboxReader {
  public:
    /**
     * Create the shoebox reader with the filename and block.
     * @param filename The file to read from
     * @param blocks The list of blocks to read
     */
    ShoeboxReader(const std::string &filename, const af::const_ref<std::size_t> &blocks)
        : filename_(filename), blocks_(blocks.begin(), blocks.end()) {
      // Check we have atleast 1 block and that block array is valid
      DIALS_ASSERT(blocks.size() >= 2);
      for (std::size_t i = 1; i < blocks.size(); ++i) {
        DIALS_ASSERT(blocks[i] > blocks[i - 1]);
      }
    }

    /**
     * @returns The filename of the shoebox file.
     */
    std::string filename() const {
      return filename_;
    }

    /**
     * @returns A list of blocks
     */
    af::shared<std::size_t> blocks() const {
      return af::shared<std::size_t>(blocks_.begin(), blocks_.end());
    }

    /**
     * @returns The number of blocks
     */
    std::size_t size() const {
      return blocks_.size() - 1;
    }

    /**
     * Return the specific block
     * @param index The index of the block
     * @returns The block z range
     */
    int2 block(std::size_t index) const {
      DIALS_ASSERT(index < size());
      return int2(blocks_[index], blocks_[index + 1]);
    }

    /**
     * Read the shoeboxes in a block
     * @param index The block index
     * @returns The list of shoeboxes in the block
     */
    af::shared<Shoebox<double> > operator[](std::size_t index) const {
      DIALS_ASSERT(index < size());
      return af::shared<Shoebox<double> >();
    }

  private:
    std::string filename_;
    af::shared<std::size_t> blocks_;
  };

  // class ShoeboxFile {
  // public:

  // ShoeboxFile(const std::string filename)
  //: filename_(filename) {

  //}

  // std::string filename() const {
  // return filename_;
  //}

  // int* block(std::size_t n) {

  //}

  // private:

  // std::string filename_;
  //};

  // class ShoeboxWriter {
  // public:

  // ShoeboxWriter(const std::string &filename
  // const af::const_ref<std::size_t> &panel,
  // const af::const_ref<int6> &bbox,
  // std::size_t max_block_size) {

  //}

  //~ShoeboxWriter() {
  // defragment();
  //}

  // void add_image(std::size_t z, std::size_t panel, const af::const_ref &image) {

  //// Ensure we can't go backward and panels are valid
  // DIALS_ASSERT(z >= current_frame_);
  // DIALS_ASSERT(panel < num_panels_);
  // current_frame_ = z;

  //// Get the current

  //}

  // private:

  // void defragment() {

  //}

  // std::size_t current_frame_;
  // std::size_t num_panels_;
  //};

}}  // namespace dials::model

#endif  // DIALS_MODEL_SERIALIZE_SHOEBOX_H
