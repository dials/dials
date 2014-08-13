

#ifndef DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H
#define DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/image.h>

namespace dials { namespace model {

  class SimpleShoeboxExtractor {
  public:

    SimpleShoeboxExtractor(const af::const_ref< Shoebox<> > shoeboxes)
      : shoeboxes_(shoeboxes.begin(), shoeboxes.end()) {

    }


    void next(const Image &image) {

    }

  private:

    af::shared< Shoebox<> > shoeboxes_;
  };

}}

#endif // DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H
