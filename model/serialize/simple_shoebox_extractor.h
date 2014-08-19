

#ifndef DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H
#define DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/image.h>

namespace dials { namespace model {

  /**
   * A class to extract shoeboxes from a sequence of images.
   */
  class SimpleShoeboxExtractor {
  public:

    struct sort_by_frame_and_panel {
      af::const_ref<int> f;
      af::const_ref<std::size_t> p;
      sort_by_frame_and_panel(
          const af::const_ref<int> frame,
          const af::const_ref<std::size_t> panel)
        : f(frame), p(panel) {}
      bool operator()(std::size_t a, std::size_t b) const {
        return (f[a] == f[b] ? p[a] < p[b] : f[a] < f[b]);
      }
    };

    /**
     * Initialise the extractor.
     * @param shoeboxes The shoeboxes to extract
     * @param frame0 The first frame
     * @param frame1 The last frame
     * @param numpanels The number of panels
     */
    SimpleShoeboxExtractor(
          const af::const_ref< Shoebox<> > shoeboxes,
          int frame0,
          int frame1,
          std::size_t numpanels)
      : shoeboxes_(shoeboxes.begin(), shoeboxes.end()),
        frame0_(frame0),
        frame1_(frame1),
        frame_(frame0),
        numframes_(frame1 - frame0),
        numpanels_(numpanels) {
      DIALS_ASSERT(frame1_ > frame0_);
      DIALS_ASSERT(numpanels > 0);

      // Compute the number of indices
      std::size_t numpartial = 0;
      for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
        int z0 = shoeboxes[i].bbox[4];
        int z1 = shoeboxes[i].bbox[5];
        DIALS_ASSERT(z1 > z0);
        numpartial += (z1 - z0);
      }

      // Check number of panels and frames
      for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
        std::size_t p = shoeboxes[i].panel;
        int z0 = shoeboxes[i].bbox[4];
        int z1 = shoeboxes[i].bbox[5];
        DIALS_ASSERT(p < numpanels_);
        DIALS_ASSERT(z0 >= frame0);
        DIALS_ASSERT(z1 <= frame1);
      }

      // Create a set of partials
      af::shared<int> frame(numpartial);
      af::shared<std::size_t> panel(numpartial);
      af::shared<std::size_t> partind(numpartial);
      af::shared<std::size_t> parent(numpartial);
      indices_.resize(numpartial);
      std::size_t j = 0;
      for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
        for (int z = shoeboxes[i].bbox[4]; z < shoeboxes[i].bbox[5]; ++z) {
          DIALS_ASSERT(j < indices_.size());
          panel[j] = shoeboxes[i].panel;
          frame[j] = z;
          partind[j] = j;
          parent[j] = i;
          j++;
        }
      }
      DIALS_ASSERT(j == numpartial);

      // Sort the indices by frame and panel
      std::sort(partind.begin(), partind.end(),
          sort_by_frame_and_panel(
            frame.const_ref(),
            panel.const_ref()));

      // Reorder indices
      indices_.resize(numpartial);
      for (std::size_t i = 0; i < partind.size(); ++i) {
        indices_[i] = parent[partind[i]];
      }

      // Create the offsets
      offset_.resize(numframes_*numpanels_+1);
      offset_[0] = 0;
      int f0 = frame0_;
      std::size_t p0 = 0;
      std::size_t k = 1;
      for (std::size_t i = 0; i < partind.size(); ++i) {
        std::size_t ii = partind[i];
        std::size_t p1 = panel[ii];
        int f1 = frame[ii];
        DIALS_ASSERT(f1 >= f0 && (f1 != f0 || p1 >= p0));
        for (; f0 < f1; ++f0) {
          for (; p0 < numpanels_; ++p0) {
            offset_[k++] = i;
          }
          p0 = 0;
        }
        for (; p0 < p1; ++p0) {
          offset_[k++] = i;
        }
        DIALS_ASSERT(f0 == f1 && p0 == p1);
      }
      for (; f0 < frame1_; ++f0) {
        for (; p0 < numpanels_; ++p0) {
          offset_[k++] = indices_.size();
        }
        p0 = 0;
      }
      DIALS_ASSERT(k == offset_.size());

      // Check the computed offsets
      for (std::size_t i = 1; i < offset_.size(); ++i) {
        DIALS_ASSERT(offset_[i] >= offset_[i-1]);
      }
      DIALS_ASSERT(offset_.back() == indices_.size());
      for (std::size_t f = 0; f < numframes_; ++f) {
        for (std::size_t p = 0; p < numpanels_; ++p) {
          int f1 = (int)f + frame0_;
          std::size_t i = p + f *numpanels_;
          std::size_t i1 = offset_[i];
          std::size_t i2 = offset_[i+1];
          DIALS_ASSERT(i2 >= i1);
          for (std::size_t i = i1; i < i2; ++i) {
            std::size_t j = indices_[i];
            std::size_t p2 = shoeboxes_[j].panel;
            int6 b = shoeboxes_[j].bbox;
            DIALS_ASSERT(f1 >= b[4] && f1 < b[5]);
            DIALS_ASSERT(p == p2);
          }
        }
      }
    }

    /**
     * Apply the next image.
     * @param image The image to apply
     */
    void next(const Image &image) {
      DIALS_ASSERT(frame_ >= frame0_ && frame_ < frame1_);
      DIALS_ASSERT(image.npanels() == numpanels_);
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref< int, af::c_grid<2> > data = image.data(p);
        af::const_ref< bool, af::c_grid<2> > mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          Shoebox<>& sbox = shoeboxes_[ind[i]];
          int6 b = sbox.bbox;
          DIALS_ASSERT(b[1] > b[0]);
          DIALS_ASSERT(b[3] > b[2]);
          DIALS_ASSERT(b[5] > b[4]);
          DIALS_ASSERT(sbox.is_consistent());
          DIALS_ASSERT(frame_ >= b[4] && frame_ < b[5]);
          int x0 = b[0];
          int x1 = b[1];
          int y0 = b[2];
          int y1 = b[3];
          int z0 = b[4];
          std::size_t xs = x1 - x0;
          std::size_t ys = y1 - y0;
          std::size_t z = frame_ - z0;
          DIALS_ASSERT(x0 >= 0 && y0 >= 0);
          DIALS_ASSERT(y1 <= data.accessor()[0]);
          DIALS_ASSERT(x1 <= data.accessor()[1]);
          for (std::size_t y = 0; y < ys; ++y) {
            for (std::size_t x = 0; x < xs; ++x) {
              sbox.data(z, y, x) = data(y+y0,x+x0);
              sbox.mask(z, y, x) = mask(y+y0,x+x0) ? Valid : 0;
            }
          }
        }
      }
      frame_++;
    }

    /**
     * @returns Is the extraction finished.
     */
    bool finished() const {
      return frame_ == frame1_;
    }

  private:

    /**
     * Get an array of indices for shoeboxes on a particular frame and panel.
     * @param frame The frame number
     * @param panel The panel number
     * @returns The indices
     */
    af::const_ref<std::size_t> indices(int frame, std::size_t panel) const {
      if (panel >= numpanels_ || frame < frame0_ || frame >= frame1_) {
        return af::const_ref<std::size_t>(0, 0);
      }
      std::size_t j0 = panel+(frame-frame0_)*numpanels_;
      DIALS_ASSERT(j0 < offset_.size()-1);
      std::size_t i0 = offset_[j0];
      std::size_t i1 = offset_[j0+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - off;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t>(&indices_[off], num);
    }

    af::shared< Shoebox<> > shoeboxes_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t numframes_;
    std::size_t numpanels_;
    af::shared<std::size_t> indices_;
    af::shared<std::size_t> offset_;
  };

}}

#endif // DIALS_MODEL_SERIALIZE_SIMPLE_SHOEBOX_EXTRACTOR_H
