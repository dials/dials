
#ifndef DIALS_PARTIAL_WRITER_H
#define DIALS_PARTIAL_WRITER_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <H5Cpp.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials {

  class PartialWriter {
  public:

    PartialWriter(
        const char *filename,
        const char *external,
        int frame0,
        int frame1,
        const af::const_ref<int> &frame,
        const af::const_ref<int4> &bbox) {
      init(filename, external, frame0, frame1, frame, bbox);
    }

    ~PartialWriter() {
      file_.close();
      std::remove(filename_);
    }

    void write(int frame, const af::const_ref< int, af::c_grid<2> > &image) {
      DIALS_ASSERT(frame >= frame0_ && frame < frame1_);
      std::size_t index = frame - frame0_;
      write_image_to_buffer(index, image);
      write_buffer_to_file(index);
    }

  private:

    H5::H5File file_;
    int frame0_;
    int frame1_;
    H5::DataSet shoeboxes_;
    std::vector<hsize_t> image_offset_;
    std::vector<std::size_t> partial_offset_;
    std::vector<int> buffer_;
    std::vector<int> frame_;
    std::vector<int4> bbox2d_;
    std::vector<std::size_t> indices_;
  };

}

#endif // DIALS_PARTIAL_WRITER_H
