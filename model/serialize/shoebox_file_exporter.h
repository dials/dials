/*
 * shoebox_file_exporter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_EXPORTER_H
#define DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_EXPORTER_H

#include <vector>
#include <stack>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <dials/model/serialize/shoebox_file.h>

namespace dials { namespace model { namespace serialize {

  using af::int2;
  using af::int6;

  /**
   * A class to buffer the shoebox data. This class finds the maximum amount of
   * memory needed at any one time, it then allocated the memory to write the
   * shoeboxes. Looping though the frames, as shoeboxes comes into scope an
   * offset into the buffer is given to the shoebox and relinqushed when the
   * shoebox goes out of scope.
   */
  class ShoeboxFileExporterBuffer {
  public:

    typedef std::vector<std::size_t> index_array_type;
    typedef std::vector<index_array_type> pf_index_array_type;

    struct sort_by_lower_z {
      af::const_ref<int6> bbox;
      sort_by_lower_z(const af::const_ref<int6> &bbox_)
        : bbox(bbox_) {}
      template <typename T>
      bool operator()(T a, T b) {
        return bbox[a][4] < bbox[b][4];
      }
    };

    struct sort_by_upper_z {
      af::const_ref<int6> bbox;
      sort_by_upper_z(const af::const_ref<int6> &bbox_)
        : bbox(bbox_) {}
      template <typename T>
      bool operator()(T a, T b) {
        return bbox[a][5] < bbox[b][5];
      }
    };

    /**
     * Initialize the buffer
     * @param bbox The bounding boxes
     * @param nframes The number of frames in the imageset
     */
    ShoeboxFileExporterBuffer(const af::const_ref<int6> &bbox, std::size_t nframes)
        : max_size_(0),
          max_num_(0),
          offset_(bbox.size()) {
      // Compute maximum size of shoebox
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const int6 &b = bbox[i];
        std::size_t size = (b[5] - b[4])*(b[3] - b[2])*(b[1] - b[0]);
        if (size > max_size_) {
          max_size_ = size;
        }
      }

      // Index of shoebox sorted by min and max z
      std::vector<std::size_t> minz(bbox.size());
      std::vector<std::size_t> maxz(bbox.size());
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        minz[i] = i;
        maxz[i] = i;
      }
      std::sort(minz.begin(), minz.end(), sort_by_lower_z(bbox));
      std::sort(maxz.begin(), maxz.end(), sort_by_upper_z(bbox));

      // Find a position available in the buffer and assign to shoebox. If no
      // position is available add another. When positions become available add
      // them back to the stack of available positions
      std::stack<std::size_t> available;
      std::size_t minz_index = 0;
      std::size_t maxz_index = 0;
      int current_frame = bbox[minz[minz_index]][4];
      while(current_frame <= nframes) {
        while (maxz_index < maxz.size()
            && bbox[maxz[maxz_index]][5] == current_frame) {
          available.push(offset_[maxz[maxz_index]]);
          maxz_index++;
        }
        while (minz_index < minz.size()
            && bbox[minz[minz_index]][4] == current_frame) {
          std::size_t pos = 0;
          if (available.empty()) {
            pos = max_num_++;
          } else {
            pos = available.top();
            available.pop();
          }
          offset_[minz[minz_index]] = pos;
          minz_index++;
        }
        current_frame++;
      }

      // Allocate the buffer
      buffer_.resize(max_num_ * max_size_, 0);
      check_.resize(max_num_, -1);
    }

    /**
     * Get the buffer for the shoebox
     */
    const int* operator[](std::size_t index) const {
      DIALS_ASSERT(index < offset_.size());
      std::size_t oindex = offset_[index];
      DIALS_ASSERT(oindex < max_num_);
      DIALS_ASSERT(check_[oindex] == index);
      return &buffer_[oindex * max_size_];
    }

    /**
     * Get the buffer for the shoebox
     */
    int* operator[](std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      std::size_t oindex = offset_[index];
      DIALS_ASSERT(oindex < max_num_);
      DIALS_ASSERT(check_[oindex] == index);
      return &buffer_[oindex * max_size_];
    }

    /**
     * Hold this chunk
     */
    void hold(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      std::size_t oindex = offset_[index];
      DIALS_ASSERT(oindex < max_num_);
      DIALS_ASSERT(check_[oindex] == -1);
      check_[oindex] = index;
    }

    /**
     * Release this chunk
     */
    void release(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      std::size_t oindex = offset_[index];
      DIALS_ASSERT(oindex < max_num_);
      DIALS_ASSERT(check_[oindex] == index);
      check_[oindex] = -1;
    }

  private:

    std::size_t max_size_;
    std::size_t max_num_;
    std::vector<int> buffer_;
    std::vector<std::size_t> offset_;
    std::vector<int> check_;
  };


  /**
   * A class to export shoeboxes to an intermediate file after extracting the
   * raw data from the image pixels.
   */
  class ShoeboxFileExporter {
  public:

    // Buffer size of 100MB
    static const std::size_t BUFFER_MAX_SIZE = 100000000 / sizeof(int);

    typedef std::vector<std::size_t> index_array_type;
    typedef std::vector<index_array_type> pf_index_array_type;
    typedef int* shoebox_type;
    typedef af::ref< int, af::c_grid<3> > shoebox_ref_type;

    /**
     * Initialise the exporter
     * @param filename The intermediate filename
     * @param panel The list of panels
     * @param bbox The list of bounding boxes
     * @param num_frame The number of frames
     * @param num_panel The number of panels
     */
    ShoeboxFileExporter(
            const std::string filename,
            const af::const_ref<std::size_t> &panel,
            const af::const_ref<int6> &bbox,
            const af::const_ref<double> &z,
            std::size_t num_frame,
            std::size_t num_panel,
            const std::string &blob)
        : writer_(filename, panel, bbox, z, BUFFER_MAX_SIZE, blob),
          bbox_(bbox.begin(), bbox.end()),
          panel_(panel.begin(), panel.end()),
          shoeboxes_(bbox, num_frame),
          num_frame_(num_frame),
          cur_frame_(0),
          num_panel_(num_panel),
          cur_panel_(0),
          per_frame_indices_(num_frame * num_panel) {

      // Ensure number of frames in valid
      DIALS_ASSERT(num_frame_ > 0);
      DIALS_ASSERT(num_panel_ > 0);

      // Initialise the shoebox indices
      init_indices();
    }

    /**
     * Flush the file to disk
     */
    void flush() {
      writer_.flush();
    }

    /**
     * @returns Is the extractor finished.
     */
    bool finished() const {
      return cur_frame_ == num_frame_;
    }

    /**
     * Add the next image. Images should be added in panel/frame order as
     * follows:
     * F0 P0
     * F0 P1
     * ...
     * F0 PN
     * F1 P0
     * F1 P1
     * ...
     * F1 PN
     * ...
     *
     * The calling code should check that the returned frame/panel pair to
     * ensure that the image data has been added as expected.
     *
     * @param image The image data
     * @returns The frame/panel numbers.
     */
    std::pair<std::size_t, std::size_t>
    next(const af::const_ref< int, af::c_grid<2> > &image) {

      // Check we're within frame and panel range
      DIALS_ASSERT(cur_frame_ < num_frame_);
      DIALS_ASSERT(cur_panel_ < num_panel_);

      // Get the list of indices of shoeboxes recorded on a frame
      const index_array_type &indices = indices_on_frame(
          cur_panel_, cur_frame_);

      // For each shoebox, add the image pixels to the shoebox, then check if
      // the shoebox has been finished. If so, then write the shoebox to disk
      // and discard the memory allocated for the shoebox
      for (std::size_t i = 0; i < indices.size(); ++i) {

        // Get some bits and pieces
        std::size_t index = indices[i];
        int6 &bbox = bbox_[index];
        std::size_t zs = bbox[5] - bbox[4];
        std::size_t ys = bbox[3] - bbox[2];
        std::size_t xs = bbox[1] - bbox[0];
        if (bbox[4] == cur_frame_) {
          shoeboxes_.hold(index);
        }
        shoebox_type shoebox = shoeboxes_[index];

        // Add the data to the shoebox
        shoebox_ref_type shoebox_ref(shoebox, af::c_grid<3>(zs, ys, xs));
        add(shoebox_ref, bbox, cur_frame_, image);

        // Write and release shoebox
        if (bbox[5] == cur_frame_+1) {
          writer_.write(index, shoebox_ref);
          shoeboxes_.release(index);
        }
      }

      // Construct the result
      std::pair<std::size_t, std::size_t> result;
      result.first = cur_frame_;
      result.second = cur_panel_;

      // Increment the panel and frame
      cur_panel_++;
      if (cur_panel_ >= num_panel_) {
        cur_panel_ = 0;
        cur_frame_++;
      }

      // Return the frame/panel pair
      return result;
    }

  private:

    /**
     * Initialise the lookup table of shoeboxes appearing on which frame/panel.
     */
    void init_indices() {
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        DIALS_ASSERT(panel_[i] < num_panel_);
        DIALS_ASSERT(bbox_[i][4] < bbox_[i][5]);
        DIALS_ASSERT(bbox_[i][4] >= 0);
        DIALS_ASSERT(bbox_[i][5] <= num_frame_);
        std::size_t p = panel_[i];
        for (std::size_t j = bbox_[i][4]; j < bbox_[i][5]; ++j) {
          per_frame_indices_[j*num_panel_ + p].push_back(i);
        }
      }
    }

    /**
     * Get the indices appearing on the panel/frame
     */
    const index_array_type& indices_on_frame(
        std::size_t panel, std::size_t frame) {
      return per_frame_indices_[frame*num_panel_ + panel];
    }

    /**
     * Add the image data to a showbox.
     * @param shoebox The shoebox to add to.
     * @param bbox The bounding box
     * @param z The frame number
     * @param image The image data
     */
    void add(shoebox_ref_type shoebox,
             const int6 &bbox,
             std::size_t z,
             const af::const_ref< int, af::c_grid<2> > &image) {

      // Check z is ok
      DIALS_ASSERT(bbox[0] >= 0);
      DIALS_ASSERT(bbox[1] <= image.accessor()[1]);
      DIALS_ASSERT(bbox[2] >= 0);
      DIALS_ASSERT(bbox[3] <= image.accessor()[0]);
      DIALS_ASSERT(z >= bbox[4] && z < bbox[5]);

      // The offsets into the image
      std::size_t k = z - bbox[4];
      std::size_t is = bbox[1] - bbox[0];
      std::size_t js = bbox[3] - bbox[2];
      std::size_t ks = bbox[5] - bbox[4];
      std::size_t i0 = bbox[0];
      std::size_t j0 = bbox[2];

      // Check shoebox is the right size
      DIALS_ASSERT(shoebox.accessor()[0] == ks);
      DIALS_ASSERT(shoebox.accessor()[1] == js);
      DIALS_ASSERT(shoebox.accessor()[2] == is);

      // Copy the image pixels
      for (std::size_t j = 0; j < js; ++j) {
        for (std::size_t i = 0; i < is; ++i) {
          shoebox(k, j, i) = image(j + j0, i + i0);
        }
      }
    }

    ShoeboxFileWriter writer_;
    af::shared<int6> bbox_;
    af::shared<std::size_t> panel_;
    ShoeboxFileExporterBuffer shoeboxes_;
    std::size_t num_frame_;
    std::size_t cur_frame_;
    std::size_t num_panel_;
    std::size_t cur_panel_;
    pf_index_array_type per_frame_indices_;
  };

}}} // namespace dials::model::serialize

#endif // DIALS_MODEL_SHOEBOX_FILE_EXPORTER_H
