/*
 * shoebox_file.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_H
#define DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_H

#include <fstream>
#include <iostream>
#include <boost/cstdint.hpp>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace model { namespace serialize {

  using boost::int32_t;
  using boost::int64_t;
  using boost::uint32_t;
  using boost::uint64_t;

  using af::int6;


  /**
   * Base class for the shoebox file that defines some flags. The shoebox file
   * is designed for efficiently reading/writing only the shoebox data to/from a
   * scratch file during processing. It is quite basic but has some structure as
   * shown below. Data is written in binary format.
   *
   * MAGIC
   * VERSION
   * HEADER_BEG
   *   N_REFLECTIONS
   *   PANEL[0]
   *   PANEL[1]
   *   ...
   *   PANEL[N-1]
   *   BBOX[0]
   *   BBOX[1]
   *   ...
   *   BBOX[N-1]
   * HEADER_END
   * DATA_BEG
   *   SHOEBOX_BEG
   *     SBOX[0][0]
   *     SBOX[0][1]
   *     ...
   *     SBOX[0][M0-1]
   *   SHOEBOX_BEG
   *     SBOX[1][0]
   *     SBOX[1][1]
   *     ...
   *     SBOX[1][M1-1]
   *   ...
   *   SHOEBOX_BEG
   *     SBOX[N-1][0]
   *     SBOX[N-1][1]
   *     ...
   *     SBOX[N-1][MNm1-1]
   * DATA_END 
   *
   * Shoeboxes are written to the file in the order of their indices, as such it
   * is up to the calling code to ensure suffient locality between shoeboxes
   * being written and subsequently read. Too much random access in the file
   * could reduce performance. It is anticipated that the bounding boxes be
   * ordered by the reflection centroid to ensure best performance when reading
   * reflections in frame order during processing. In this scheme there will be
   * some seeking during writing as shoeboxes will generally be extracted from
   * the image data out of order.
   */
  class ShoeboxFileBase {
  public:
    const static uint32_t MAGIC = 58008;
    const static uint32_t VERSION = 1;
    const static uint32_t HEADER_BEG = 2;
    const static uint32_t HEADER_END = 3;
    const static uint32_t DATA_BEG = 4;
    const static uint32_t DATA_END = 5;
    const static uint32_t SHOEBOX_BEG = 6;
  };


  /**
   * A class to efficiently write shoeboxes to file.
   */
  class ShoeboxFileWriter : public ShoeboxFileBase {
  public:

    /**
     * Setup the shoebox file. This function will initialize the file to the
     * correct size and will write the header information. Shoeboxes will be
     * written to file in the order in which the bounding boxes have been
     * specified.
     * @param filename The file to write to
     * @param bbox The list of bounding boxes.
     */
    ShoeboxFileWriter(const std::string &filename, 
                      const af::const_ref<std::size_t> &panel,
                      const af::const_ref<int6> &bbox)
        : panel_(panel.begin(), panel.end()),
          bbox_(bbox.begin(), bbox.end()),
          offset_(bbox_.size() + 1),
          file_(filename.c_str(), 
                std::ios_base::binary | std::ios_base::trunc) {
      DIALS_ASSERT(panel.size() == bbox.size());
      write_internal(MAGIC);
      write_internal(VERSION);
      write_header();
      write_data();
    }

    /**
     * Flush the data to disk.
     */
    void flush() {
      file_.flush();
      DIALS_ASSERT(file_.good());
    }

    /**
     * Write a shoebox to the file.
     * @param index The index of the shoebox
     * @param data The shoebox data
     */
    void write(
        std::size_t index, 
        const af::const_ref< int, af::c_grid<3> > &data) {
      
      // Ensure the index is in range
      DIALS_ASSERT(index < bbox_.size());

      // Ensure the size matches the bounding box
      std::size_t zs = bbox_[index][5] - bbox_[index][4];
      std::size_t ys = bbox_[index][3] - bbox_[index][2];
      std::size_t xs = bbox_[index][1] - bbox_[index][0];
      DIALS_ASSERT(data.accessor()[0] == zs);
      DIALS_ASSERT(data.accessor()[1] == ys);
      DIALS_ASSERT(data.accessor()[2] == xs);

      // Move to the desired position in the file
      uint64_t offset = data_offset_ + offset_[index];
      seek_internal(offset, std::ios_base::beg);

      // Write the data and the shoebox flag
      write_internal(SHOEBOX_BEG);
      write_internal(&data[0], data.size());
    }

  private:

    /**
     * Write the header information
     */
    void write_header() {
      
      // Begin the header
      write_internal(HEADER_BEG);
    
      // Write the panels and bboxes
      write_internal((uint32_t)bbox_.size());
      for (std::size_t i = 0; i < panel_.size(); ++i) {
        write_internal((uint32_t)panel_[i]);
      }
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
          write_internal((int32_t)bbox_[i][j]);
        }
      }

      // Calculate the offsets
      offset_[0] = 0;
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        uint64_t box_size = 
          (bbox_[i][1] - bbox_[i][0]) *
          (bbox_[i][3] - bbox_[i][2]) *
          (bbox_[i][5] - bbox_[i][4]);
        offset_[i+1] = offset_[i] + box_size * sizeof(int) + sizeof(SHOEBOX_BEG); 
      }

      // End the header
      write_internal(HEADER_END);
    }

    /**
     * Write the shoebox data
     */
    void write_data() {
      
      // Allocate space for the data
      write_internal(DATA_BEG);
      data_offset_ = file_.tellp();
      seek_internal(offset_.back(), std::ios_base::cur);
      write_internal(DATA_END);
    }

    /**
     * Write the value and check the file flags
     */
    template <typename T>
    void write_internal(T v) {
      file_.write((const char *)&v, sizeof(T));
      DIALS_ASSERT(file_.good());
    }

    /**
     * Write the block of data and check the file flags
     */
    template <typename T>
    void write_internal(const T *v, std::size_t num) {
      file_.write((const char *)v, num * sizeof(T));
      DIALS_ASSERT(file_.good());
    }

    /**
     * Seek to the desired point and check the file flags
     */
    void seek_internal(std::streamoff off, std::ios_base::seekdir way) {
      file_.seekp(off, way);
      DIALS_ASSERT(file_.good());
    }

    af::shared<int6> bbox_;
    af::shared<std::size_t> panel_;
    af::shared<uint64_t> offset_;
    std::ofstream file_;
    uint64_t data_offset_;
  };


  /**
   * A class to read the shoebox file.
   */
  class ShoeboxFileReader : public ShoeboxFileBase {
  public:
    
    /**
     * Read in the shoebox header info.
     * @param filename The file to read
     */
    ShoeboxFileReader(const std::string filename) 
        : file_(filename.c_str()) {
    
      // Read and check magic and version numbers
      DIALS_ASSERT(read_internal<uint32_t>() == MAGIC);
      DIALS_ASSERT(read_internal<uint32_t>() == VERSION);

      // Read and check the header
      check_header();

      // Read and check the data
      check_data();
    }

    /**
     * @returns The bounding boxes
     */
    af::shared<int6> bbox() const {
      return af::shared<int6>(bbox_.begin(), bbox_.end());
    }

    /**
     * @returns The panels
     */
    af::shared<std::size_t> panels() const {
      return af::shared<std::size_t>(panel_.begin(), panel_.end());
    }

    /**
     * Read a shoebox
     * @param index The index of the shoebox to read
     * @returns The shoebox data.
     */
    af::versa< int, af::c_grid<3> > read(std::size_t index) {

      // Check the index is valid.
      DIALS_ASSERT(index < bbox_.size());

      // Allocate the memory for the shoebox
      std::size_t zs = bbox_[index][5] - bbox_[index][4];
      std::size_t ys = bbox_[index][3] - bbox_[index][2];
      std::size_t xs = bbox_[index][1] - bbox_[index][0];
      af::versa< int, af::c_grid<3> > data(
          af::c_grid<3>(zs, ys, xs));

      // Move to the desired position
      uint64_t offset = data_offset_ + offset_[index];
      seek_internal(offset, std::ios_base::beg);

      // Read the data
      DIALS_ASSERT(read_internal<uint32_t>() == SHOEBOX_BEG);
      read_internal(&data[0], data.size());

      // Return the array
      return data;
    }

  private:

    /**
     * Read the header information.
     */
    void check_header() {
    
      // Check the header begin flag
      DIALS_ASSERT(read_internal<uint32_t>() == HEADER_BEG);

      // Read the number of bboxes
      uint32_t bbox_size = read_internal<uint32_t>();
      panel_.resize(bbox_size);
      bbox_.resize(bbox_size);
      offset_.resize(bbox_size + 1); 
      for (std::size_t i = 0; i < bbox_size; ++i) {
        panel_[i] = read_internal<uint32_t>();
      }
      for (std::size_t i = 0; i < bbox_size; ++i) {
        int6 b;
        for (std::size_t j = 0; j < 6; ++j) {
          b[j] = read_internal<int32_t>(); 
        }
        DIALS_ASSERT(b[1] > b[0]);
        DIALS_ASSERT(b[3] > b[2]);
        DIALS_ASSERT(b[5] > b[4]);
        bbox_[i] = b;
      }

      // Calculate the offsets
      offset_[0] = 0;
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        uint64_t box_size = 
          (bbox_[i][1] - bbox_[i][0]) *
          (bbox_[i][3] - bbox_[i][2]) *
          (bbox_[i][5] - bbox_[i][4]);
        offset_[i+1] = offset_[i] + box_size * sizeof(int) + sizeof(SHOEBOX_BEG); 
      }

      // Check the header end flag
      DIALS_ASSERT(read_internal<uint32_t>() == HEADER_END);
    }

    /**
     * Read the data and check the shoebox flags are correct.
     */
    void check_data() {

      // Check the data flags
      DIALS_ASSERT(read_internal<uint32_t>() == DATA_BEG);
      data_offset_ = file_.tellg();

      // Check all the shoebox data flags
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        DIALS_ASSERT(read_internal<uint32_t>() == SHOEBOX_BEG);
        seek_internal(data_offset_ + offset_[i+1], std::ios_base::beg);
      }

      //file_.seekg(offset_.back(), std::ios_base::cur);
      DIALS_ASSERT(read_internal<uint32_t>() == DATA_END);
     
      // Try to read a byte and check eof
      read_internal_nocheck<char>();
      DIALS_ASSERT(file_.eof());

      // Clear EOF and go to start of data
      file_.clear();
      seek_internal(data_offset_, std::ios_base::beg);
    }

    /**
     * Read the value and check the file flags
     */
    template <typename T>
    T read_internal() {
      T value;
      file_.read((char *)&value, sizeof(T));
      DIALS_ASSERT(file_.good());
      return value;
    }

    /**
     * Read the value and check the file flags
     */
    template <typename T>
    T read_internal_nocheck() {
      T value;
      file_.read((char *)&value, sizeof(T));
      return value;
    }

    /**
     * Read the data block and check the flags.
     */
    template <typename T>
    void read_internal(T *data, std::size_t num) {
      file_.read((char *)data, num * sizeof(T));
      DIALS_ASSERT(file_.good());
    }

    /**
     * Seek and check the file flags.
     */
    void seek_internal(std::streamoff off, std::ios_base::seekdir way) {
      file_.seekg(off, way);
      DIALS_ASSERT(file_.good());
    }

    std::ifstream file_;
    af::shared<int6> bbox_;
    af::shared<std::size_t> panel_;
    af::shared<uint64_t> offset_;
    uint64_t data_offset_;
  };

}}} // namespace dials::model::serialize

#endif // DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_H
