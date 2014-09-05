/*
 * interface.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
#define DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <boost/function.hpp>
#include <cctbx/miller/index_generator.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>


namespace dials { namespace algorithms {

  using model::Image;
  using model::Shoebox;

  namespace detail {

    /**
     * A class to extract shoebox pixels from images
     */
    class ShoeboxExtractor {
    public:

      ShoeboxExtractor() {}

      /**
       * Initialise the index array. Determine which reflections are recorded on
       * each frame and panel ahead of time to enable quick lookup of the
       * reflections to be written to when processing each image.
       */
      ShoeboxExtractor(af::reflection_table data,
                       std::size_t npanels,
                       int frame0,
                       int frame1)
          : npanels_(npanels),
            frame0_(frame0),
            frame1_(frame1) {
        DIALS_ASSERT(frame0_ < frame1_);
        DIALS_ASSERT(npanels_ > 0);
        DIALS_ASSERT(data.is_consistent());
        DIALS_ASSERT(data.contains("panel"));
        DIALS_ASSERT(data.contains("bbox"));
        DIALS_ASSERT(data.contains("shoebox"));
        shoebox_ = data["shoebox"];
        std::size_t nframes = frame1_ - frame0_;
        af::const_ref<std::size_t> panel = data["panel"];
        af::const_ref<int6> bbox = data["bbox"];
        std::size_t size = nframes * npanels_;
        std::vector<std::size_t> num(size, 0);
        std::vector<std::size_t> count(size, 0);
        for (std::size_t i = 0; i < bbox.size(); ++i) {
          DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
          DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
          DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
          for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
            std::size_t j = panel[i] + (z - frame0_)*npanels_;
            DIALS_ASSERT(j < num.size());
            num[j]++;
          }
        }
        offset_.push_back(0);
        std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));
        indices_.resize(offset_.back());
        for (std::size_t i = 0; i < bbox.size(); ++i) {
          for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
            std::size_t j = panel[i] + (z - frame0_)*npanels_;
            std::size_t k = offset_[j] + count[j];
            DIALS_ASSERT(j < count.size());
            DIALS_ASSERT(k < indices_.size());
            indices_[k] = i;
            count[j]++;
          }
        }
        DIALS_ASSERT(count == num);
      }

      /**
       * Extract the pixels from the image and copy to the relevant shoeboxes.
       * @param image The image to process
       * @param frame The current image frame
       */
      void next(const Image &image, int frame) {
        typedef Shoebox<>::float_type float_type;
        typedef af::ref<float_type, af::c_grid<3> > sbox_data_type;
        typedef af::ref<int,        af::c_grid<3> > sbox_mask_type;
        DIALS_ASSERT(frame >= frame0_ && frame < frame1_);
        DIALS_ASSERT(image.npanels() == npanels_);
        for (std::size_t p = 0; p < image.npanels(); ++p) {
          af::const_ref<std::size_t> ind = indices(frame, p);
          af::const_ref< int, af::c_grid<2> > data = image.data(p);
          af::const_ref< bool, af::c_grid<2> > mask = image.mask(p);
          DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
          for (std::size_t i = 0; i < ind.size(); ++i) {
            DIALS_ASSERT(ind[i] < shoebox_.size());
            Shoebox<>& sbox = shoebox_[ind[i]];
            int6 b = sbox.bbox;
            sbox_data_type sdata = sbox.data.ref();
            sbox_mask_type smask = sbox.mask.ref();
            DIALS_ASSERT(b[1] > b[0]);
            DIALS_ASSERT(b[3] > b[2]);
            DIALS_ASSERT(b[5] > b[4]);
            DIALS_ASSERT(frame >= b[4] && frame < b[5]);
            int x0 = b[0];
            int x1 = b[1];
            int y0 = b[2];
            int y1 = b[3];
            int z0 = b[4];
            std::size_t xs = x1 - x0;
            std::size_t ys = y1 - y0;
            std::size_t z = frame - z0;
            DIALS_ASSERT(x0 >= 0 && y0 >= 0);
            DIALS_ASSERT(y1 <= data.accessor()[0]);
            DIALS_ASSERT(x1 <= data.accessor()[1]);
            DIALS_ASSERT(sbox.is_consistent());
            for (std::size_t y = 0; y < ys; ++y) {
              for (std::size_t x = 0; x < xs; ++x) {
                sdata(z, y, x) = data(y+y0,x+x0);
                smask(z, y, x) = mask(y+y0,x+x0) ? Valid : 0;
              }
            }
          }
        }
      }

    private:

      /**
       * Get an index array specifying which reflections are recorded on a given
       * frame and panel.
       * @param frame The frame number
       * @param panel The panel number
       * @returns An array of indices
       */
      af::const_ref<std::size_t> indices(int frame, std::size_t panel) const {
        std::size_t j0 = panel+(frame-frame0_)*npanels_;
        DIALS_ASSERT(offset_.size() > 0);
        DIALS_ASSERT(j0 < offset_.size()-1);
        std::size_t i0 = offset_[j0];
        std::size_t i1 = offset_[j0+1];
        DIALS_ASSERT(i1 >= i0);
        std::size_t off = i0;
        std::size_t num = i1 - off;
        DIALS_ASSERT(off + num <= indices_.size());
        return af::const_ref<std::size_t>(&indices_[off], num);
      }

      std::size_t npanels_;
      int frame0_;
      int frame1_;
      af::shared< Shoebox<> > shoebox_;
      std::vector<std::size_t> indices_;
      std::vector<std::size_t> offset_;
    };

  }

  /**
   * A class to extract shoeboxes from a sequence of images.
   */
  class IntegrationTask3DMultiExecutor {
  public:

    /**
     * Initialise the extractor.
     * @param shoeboxes The shoeboxes to extract
     * @param frame0 The first frame
     * @param frame1 The last frame
     * @param npanels The number of panels
     */
    IntegrationTask3DMultiExecutor(
          af::reflection_table data,
          tiny<int,2> job,
          std::size_t npanels)
      : data_(data),
        frame0_(job[0]),
        frame1_(job[1]),
        frame_(frame0_),
        nframes_(frame1_ - frame0_),
        npanels_(npanels) {
      DIALS_ASSERT(frame1_ > frame0_);
      DIALS_ASSERT(npanels > 0);
      extractor_ = detail::ShoeboxExtractor(data, npanels, frame0_, frame1_);
    }

    /** @returns The first frame.  */
    int frame0() const {
      return frame0_;
    }

    /** @returns The last frame */
    int frame1() const {
      return frame1_;
    }

    /** @returns The current frame. */
    int frame() const {
      return frame_;
    }

    /** @returns The number of frames  */
    std::size_t nframes() const {
      return nframes_;
    }

    /** @returns The number of panels */
    std::size_t npanels() const {
      return npanels_;
    }

    /**
     * Apply the next image.
     * @param image The image to apply
     */
    void next(const Image &image) {
      extractor_.next(image, frame_);
      frame_++;
    }

    /**
     * @returns Is the extraction finished.
     */
    bool finished() const {
      return frame_ == frame1_;
    }

    /**
     * @returns The integrated reflection data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return data_;
    }

  private:

    af::reflection_table data_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    std::size_t npanels_;
    detail::ShoeboxExtractor extractor_;
  };

  namespace detail {

    /**
     * A class to specify the 3D integration task.
     */
    class JobManager {
    public:

      /**
       * Initialise the task
       * @param data The reflection data
       * @param npanels The number of panels
       * @param jobs The list of integration jobs
       * @param off The offset into the ind array for each job
       * @param ind The indices for each job
       */
      JobManager(
              const af::const_ref<std::size_t> &panel,
              const af::const_ref<int6> &bbox,
              const af::const_ref<std::size_t> &flags,
              const af::const_ref< tiny<int,2> > &jobs,
              std::size_t npanels)
          : jobs_(jobs.begin(), jobs.end()),
            npanels_(npanels),
            ndata_(panel.size()) {

        // Check some input
        DIALS_ASSERT(panel.size() == bbox.size());
        DIALS_ASSERT(panel.size() == flags.size());
        DIALS_ASSERT(jobs.size() > 0);
        DIALS_ASSERT(npanels > 0);

        // Compute the job indices
        compute_indices(panel, bbox, flags);

        // Check the jobs are valid. Jobs must be ordered in increasing frame
        // number and can overlap but 1 jobs must not be fully contained in
        // another. There must also be no gaps.
        for (std::size_t i = 0; i < jobs.size()-1; ++i) {
          DIALS_ASSERT(jobs[i][0] <  jobs[i][1]);
          DIALS_ASSERT(jobs[i][0] <  jobs[i+1][0]);
          DIALS_ASSERT(jobs[i][1] <  jobs[i+1][1]);
          DIALS_ASSERT(jobs[i][1] >= jobs[i+1][0]);
        }

        // Get the range of frames
        frame0_ = jobs_.front()[0];
        frame1_ = jobs_.back()[1];
        DIALS_ASSERT(frame0_ < frame1_);
        nframes_ = frame1_ - frame0_;

        // Check the jobs ids are valid
        DIALS_ASSERT(offset_.size() == njobs()+1);
        DIALS_ASSERT(offset_.front() == 0);
        DIALS_ASSERT(offset_.back() == indices_.size());
        for (std::size_t i = 0; i < njobs(); ++i) {
          af::const_ref<std::size_t> ind = indices(i);
          for (std::size_t j = 0; j < ind.size(); ++j) {
            std::size_t k = ind[j];
            DIALS_ASSERT(k < panel.size());
            std::size_t p = panel[k];
            int z0 = bbox[k][4];
            int z1 = bbox[k][5];
            DIALS_ASSERT(p < npanels_);
            DIALS_ASSERT(z1 > z0);
            int j0 = job(i)[0];
            int j1 = job(i)[1];
            DIALS_ASSERT(j1 > j0);
            DIALS_ASSERT(z0 >= j0 && z1 <= j1);
          }
        }
      }

      /**
       * @returns The job frame range
       */
      tiny<int,2> job(std::size_t index) const {
        DIALS_ASSERT(index < jobs_.size());
        return jobs_[index];
      }

      /**
       * @param index The index of the job
       * @returns The indices of reflection on a job
       */
      af::const_ref<std::size_t> indices(std::size_t index) const {
        DIALS_ASSERT(index < offset_.size()-1);
        std::size_t i0 = offset_[index];
        std::size_t i1 = offset_[index+1];
        DIALS_ASSERT(i1 >= i0);
        std::size_t off = i0;
        std::size_t num = i1 - i0;
        DIALS_ASSERT(off + num <= indices_.size());
        return af::const_ref<std::size_t> (&indices_[off], num);
      }

      /**
       * @param index The index of the job
       * @returns The mask of reflections on a job
       */
      af::const_ref<bool> mask(std::size_t index) const {
        DIALS_ASSERT(index < offset_.size()-1);
        std::size_t i0 = offset_[index];
        std::size_t i1 = offset_[index+1];
        DIALS_ASSERT(i1 >= i0);
        std::size_t off = i0;
        std::size_t num = i1 - i0;
        DIALS_ASSERT(off + num <= mask_.size());
        return af::const_ref<bool> (&mask_[off], num);
      }

      /**
       * @returns The first frame in the task.
       */
      int frame0() const {
        return frame0_;
      }

      /**
       * @returns The last frame in the task
       */
      int frame1() const {
        return frame1_;
      }

      /**
       * @returns The number of frames in the task
       */
      std::size_t nframes() const {
        return nframes_;
      }

      /**
       * @returns The number of panels
       */
      std::size_t npanels() const {
        return npanels_;
      }

      /**
       * @returns The number of jobs
       */
      std::size_t njobs() const {
        return jobs_.size();
      }

      /**
       * @returns The number of data items
       */
      std::size_t ndata() const {
        return ndata_;
      }

    private:

      /**
       * Compute the indices of the reflections in each job.
       */
      void compute_indices(
          const af::const_ref<std::size_t> &panel,
          const af::const_ref<int6> &bbox,
          const af::const_ref<std::size_t> &flags) {

        typedef std::pair<std::size_t, bool> job_type;
        typedef std::vector<job_type> job_list_type;

        // Get which reflections to process in which job
        std::vector<job_list_type> indices(jobs_.size());
        for (std::size_t index = 0; index < bbox.size(); ++index) {
          int z0 = bbox[index][4];
          int z1 = bbox[index][5];
          std::size_t f = flags[index];
          DIALS_ASSERT(!(f & af::DontIntegrate));
          std::size_t jmin = 0;
          double dmin = 0;
          bool first = true;
          for (std::size_t j = 0; j < jobs_.size(); ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              if (f & af::ReferenceSpot) {
                indices[j].push_back(job_type(index, false));
              }
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (first || d < dmin) {
                jmin = j;
                dmin = d;
                first = false;
              }
            }
          }
          DIALS_ASSERT(first != true);
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          DIALS_ASSERT(z0 >= jz0 && z1 <= jz1);
          if (f & af::ReferenceSpot) {
            DIALS_ASSERT(indices[jmin].size() > 0);
            DIALS_ASSERT(indices[jmin].back().first == index);
            indices[jmin].back().second = true;
          } else {
            indices[jmin].push_back(job_type(index, true));
          }
        }

        // Compute number of reflections in each task
        std::vector<std::size_t> num(jobs_.size(), 0);
        for (std::size_t i = 0; i < num.size(); ++i) {
          num[i] = indices[i].size();
        }

        // Compute offsests
        offset_.push_back(0);
        std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));

        // Compute indices
        indices_.resize(offset_.back());
        mask_.resize(offset_.back());
        std::size_t k = 0;
        for (std::size_t i = 0; i < indices.size(); ++i) {
          const job_list_type& ind = indices[i];
          for (std::size_t j = 0; j < ind.size(); ++j, ++k) {
            DIALS_ASSERT(k < indices_.size());
            indices_[k] = ind[j].first;
            mask_[k] = ind[j].second;
            DIALS_ASSERT(indices_[k] < bbox.size());
          }
        }
        DIALS_ASSERT(k == indices_.size());
      }

      std::vector< tiny<int,2> > jobs_;
      std::size_t npanels_;
      std::vector< std::size_t > offset_;
      std::vector< std::size_t > indices_;
      af::shared< bool > mask_;
      int frame0_;
      int frame1_;
      std::size_t nframes_;
      std::size_t ndata_;
    };


    /**
     * A class to help in allocating memory for 3D shoeboxes
     *
     * FIXME Add custom allocator for shoebox data and do a single
     * allocation at the start of the process
     */
    class ShoeboxAllocator {
    public:

      /**
       * Initialise the allocator with the task spec
       * @param spec The integration task specification
       */
      ShoeboxAllocator(const JobManager &job_manager) {

        typedef std::pair<std::size_t,std::size_t> active_type;

        // Determine the first and last job on which the reflections
        // are active
        active_type none(job_manager.njobs(), 0);
        std::vector<active_type> active(job_manager.ndata(), none);
        for (std::size_t i = 0; i < job_manager.njobs(); ++i) {
          af::const_ref<std::size_t> ind = job_manager.indices(i);
          for (std::size_t j = 0; j < ind.size(); ++j) {
            std::size_t k = ind[j];
            DIALS_ASSERT(k < active.size());
            active_type &a = active[k];
            if (i < a.first)  a.first  = i;
            if (i > a.second) a.second = i;
          }
        }

        // Some helper arrays
        std::vector<std::size_t> mcount(job_manager.njobs(), 0);
        std::vector<std::size_t> fcount(job_manager.njobs(), 0);
        std::vector<std::size_t> mnum(job_manager.njobs(), 0);
        std::vector<std::size_t> fnum(job_manager.njobs(), 0);

        // Determine the number of reflections that need to be allocated
        // and free at the beginning and end of each job
        for (std::size_t i = 0; i < active.size(); ++i) {
          DIALS_ASSERT(active[i].first < mnum.size());
          DIALS_ASSERT(active[i].second < fnum.size());
          mnum[active[i].first]++;
          fnum[active[i].second]++;
        }

        // Resize the arrays
        moffset_.resize(job_manager.njobs() + 1);
        foffset_.resize(job_manager.njobs() + 1);
        mindices_.resize(job_manager.ndata());
        findices_.resize(job_manager.ndata());

        // Compute the offsets
        moffset_[0] = 0;
        foffset_[0] = 0;
        std::partial_sum(mnum.begin(), mnum.end(), moffset_.begin()+1);
        std::partial_sum(fnum.begin(), fnum.end(), foffset_.begin()+1);

        // Check the offsets are correct
        DIALS_ASSERT(moffset_.back() == mindices_.size());
        DIALS_ASSERT(foffset_.back() == findices_.size());

        // Compute the indices
        for (std::size_t i = 0; i < active.size(); ++i) {
          std::size_t i1 = active[i].first;
          std::size_t i2 = active[i].second;
          DIALS_ASSERT(i1 < moffset_.size());
          DIALS_ASSERT(i2 < foffset_.size());
          std::size_t k1 = moffset_[i1] + mcount[i1];
          std::size_t k2 = foffset_[i2] + fcount[i2];
          DIALS_ASSERT(k1 < mindices_.size());
          DIALS_ASSERT(k2 < findices_.size());
          mindices_[k1] = i;
          findices_[k2] = i;
          mcount[i1]++;
          fcount[i2]++;
        }

        // Check the counts are as expected
        DIALS_ASSERT(mcount == mnum);
        DIALS_ASSERT(fcount == fnum);
      }

      /**
       * @param index The job index
       * @returns The indices of reflections to allocate
       */
      af::const_ref<std::size_t> malloc(std::size_t index) const {
        DIALS_ASSERT(index < moffset_.size()-1);
        std::size_t i0 = moffset_[index];
        std::size_t i1 = moffset_[index+1];
        DIALS_ASSERT(i1 >= i0);
        std::size_t off = i0;
        std::size_t num = i1 - i0;
        DIALS_ASSERT(off + num <= mindices_.size());
        return af::const_ref<std::size_t> (&mindices_[off], num);
      }

      /**
       * @param index The job index
       * @returns The indices of reflections to free
       */
      af::const_ref<std::size_t> free(std::size_t index) const {
        DIALS_ASSERT(index < foffset_.size()-1);
        std::size_t i0 = foffset_[index];
        std::size_t i1 = foffset_[index+1];
        DIALS_ASSERT(i1 >= i0);
        std::size_t off = i0;
        std::size_t num = i1 - i0;
        DIALS_ASSERT(off + num <= findices_.size());
        return af::const_ref<std::size_t> (&findices_[off], num);
      }

    private:

      std::vector<std::size_t> moffset_;
      std::vector<std::size_t> foffset_;
      std::vector<std::size_t> mindices_;
      std::vector<std::size_t> findices_;
    };

  }


  /**
   * A class to execute an integration task
   */
  class IntegrationTask3DExecutor {
  public:

    typedef af::reflection_table rtable;
    typedef boost::function<rtable (rtable)> callback_type;

    /**
     * Initialise the exector
     * @param spec The integration specification
     * @param callback The function to call to process the data
     */
    IntegrationTask3DExecutor(
            af::reflection_table data,
            const af::const_ref< tiny<int,2> > &jobs,
            std::size_t npanels,
            callback_type callback)
        : data_(data),
          job_manager_(
              data["panel"],
              data["bbox"],
              data["flags"],
              jobs,
              npanels),
          allocator_(job_manager_),
          process_(callback) {

      // Check the reflection table contains a few required fields
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("panel"));
      DIALS_ASSERT(data.contains("bbox"));

      // Set the starting frame and current active job
      frame_ = job_manager_.frame0();
      begin_active_ = 0;
      end_active_ = 0;

      // Initialise the shoeboxes
      af::const_ref< std::size_t > panel = data_["panel"];
      af::const_ref< int6 > bbox = data_["bbox"];
      shoebox_ = data_["shoebox"];
      for (std::size_t i = 0; i < shoebox_.size(); ++i) {
        shoebox_[i] = Shoebox<>(panel[i], bbox[i]);
      }

      // Initialise the offsets and indices for each frame/panel
      extractor_ = detail::ShoeboxExtractor(data_, npanels, frame0(), frame1());
    }

    /**
     * @returns The first frame.
     */
    int frame0() const {
      return job_manager_.frame0();
    }

    /**
     * @returns The last frame
     */
    int frame1() const {
      return job_manager_.frame1();
    }

    /**
     * @returns The current frame.
     */
    int frame() const {
      return frame_;
    }

    /**
     * @returns The number of frames
     */
    std::size_t nframes() const {
      return job_manager_.nframes();
    }

    /**
     * @returns The number of panels
     */
    std::size_t npanels() const {
      return job_manager_.npanels();
    }

    /**
     * Process the next image
     * @param image The image to process
     */
    void next(const Image &image) {
      DIALS_ASSERT(!finished());
      if (begin_active_ < job_manager_.njobs()) {
        int j0 = job_manager_.job(begin_active_)[0];
        DIALS_ASSERT(frame_ <= j0);
        if (frame_ == j0) {
          af::const_ref<std::size_t> ind = allocator_.malloc(begin_active_);
          for (std::size_t i = 0; i < ind.size(); ++i) {
            shoebox_[ind[i]].allocate();
          }
          begin_active_++;
        }
      }
      extractor_.next(image, frame_);
      frame_++;
      int j0 = job_manager_.job(end_active_)[0];
      int j1 = job_manager_.job(end_active_)[1];
      DIALS_ASSERT(frame_ >= j0 && frame_ <= j1);
      if (frame_ == j1) {
        merge(process_(split(end_active_)), end_active_);
        af::const_ref<std::size_t> ind = allocator_.free(end_active_);
        for (std::size_t i = 0; i < ind.size(); ++i) {
          shoebox_[ind[i]].deallocate();
        }
        end_active_++;
      }
    }

    /**
     * @returns The integrated reflection data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      data_.erase("shoebox");
      return data_;
    }

    /**
     * @returns Is the task finished
     */
    bool finished() const {
      return frame_ == frame1();
    }

  private:

    /**
     * Select the reflections for the current job
     */
    af::reflection_table split(std::size_t job) const {
      using namespace dials::af::boost_python::flex_table_suite;
      af::const_ref<bool> mask = job_manager_.mask(job);
      af::reflection_table result = select_rows_index(
          data_,
          job_manager_.indices(job));
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (!mask[i]) {
          flags[i] |= af::DontIntegrate;
        }
      }
      return result;
    }

    /**
     * Merge the results from the current job
     */
    void merge(af::reflection_table result, std::size_t job) {
      using namespace dials::af::boost_python::flex_table_suite;
      set_selected_rows_index_mask(
          data_,
          job_manager_.indices(job),
          job_manager_.mask(job),
          result);
    }

    af::reflection_table data_;
    detail::JobManager job_manager_;
    detail::ShoeboxAllocator allocator_;
    detail::ShoeboxExtractor extractor_;
    int frame_;
    std::size_t begin_active_;
    std::size_t end_active_;
    af::shared< Shoebox<> > shoebox_;
    callback_type process_;
  };

  namespace detail {

    /**
     * Compute the integration jobs
     * @param array_range The range of frames to process
     * @param block_size The number of frames in a job
     */
    af::shared< tiny<int,2> > compute_jobs(
          tiny<int,2> array_range,
          double block_size) {
      int frame0 = array_range[0];
      int frame1 = array_range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      DIALS_ASSERT(nframes >= 2);
      if (block_size > nframes) {
        block_size = nframes;
      }
      DIALS_ASSERT(block_size >= 2);
      int nblocks = (int)std::ceil(2.0 * nframes / block_size);
      DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
      int half_block_size = (int)std::ceil((double)nframes / (double)nblocks);
      af::shared<int> indices;
      indices.push_back(frame0);
      for (int i = 0; i < nblocks; ++i) {
        int frame = frame0 + (i + 1) * half_block_size;
        if (frame > frame1) {
          frame = frame1;
        }
        indices.push_back(frame);
        if (frame == frame1) {
          break;
        }
      }
      DIALS_ASSERT(indices.front() == frame0);
      DIALS_ASSERT(indices.back() == frame1);
      DIALS_ASSERT(indices.size() > 2);
      af::shared< tiny<int,2> > jobs;
      for (std::size_t i = 0; i < indices.size() - 2; ++i) {
        int i1 = indices[i];
        int i2 = indices[i+2];
        DIALS_ASSERT(i2 > i1);
        jobs.push_back(tiny<int,2>(i1, i2));
      }
      DIALS_ASSERT(jobs.size() > 0);
      return jobs;
    }

  }

  /**
   * A class to do the integration management
   */
  class IntegrationManager3DExecutor {
  public:

    /**
     * Initialise the manager executor class.
     * @param data The reflection data
     * @param array_range The array range
     * @param block_size The size of the blocks
     * @param num_tasks The number of tasks to make
     * @param npanels The number of panels
     */
    IntegrationManager3DExecutor(
        af::reflection_table data,
        vec2<int> array_range,
        double block_size)
          : data_(data),
            array_range_(array_range),
            block_size_(block_size),
            finished_(false) {
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(array_range[1] > array_range[0]);
      DIALS_ASSERT(block_size > 0);
      jobs_ = detail::compute_jobs(array_range, block_size);
      compute_indices();
    }

    /**
     * @returns The reflection data
     */
    af::reflection_table data() const {
      DIALS_ASSERT(finished());
      return data_;
    }

    /**
     * @returns The number of tasks
     */
    std::size_t size() const {
      return 1;
    }

    /**
     * @returns Have all the tasks completed.
     */
    bool finished() const {
      return finished_;
    }

    /**
     * @returns The list of jobs
     */
    af::shared< tiny<int,2> > job(std::size_t index) const {
      DIALS_ASSERT(index == 0);
      return jobs_;
    }

    /**
     * @returns A list of ignored reflections
     */
    af::shared<std::size_t> ignored() const {
      return af::shared<std::size_t>(ignored_);
    }

    /**
     * Get the specification of the requested task
     * @param index The task index
     * @returns The task spec
     */
    af::reflection_table split(std::size_t index) const {
      using namespace af::boost_python::flex_table_suite;
      DIALS_ASSERT(index == 0);
      return select_rows_index(data_, process_.const_ref());
    }

    /**
     * Accumulate the results from the separate tasks
     * @param index The index of the task
     * @param result The result of the task
     */
    void accumulate(std::size_t index, af::reflection_table result) {
      using namespace af::boost_python::flex_table_suite;
      DIALS_ASSERT(index == 0);
      DIALS_ASSERT(!finished_);
      DIALS_ASSERT(result.size() == process_.size());
      set_selected_rows_index(data_, process_.const_ref(), result);
      finished_ = true;
    }

  private:

    /**
     * Compute the indices of the lookup tables
     */
    void compute_indices() {

      typedef std::pair<std::size_t, bool> job_type;
      typedef std::vector<job_type> job_list_type;

      // Get some reflection data
      af::const_ref<int6> bbox = data_["bbox"];
      af::ref<std::size_t> flags = data_["flags"];

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t jmin = 0;
          double dmin = 0;
          bool first = true;
          for (std::size_t j = 0; j < jobs_.size(); ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (first || d < dmin) {
                jmin = j;
                dmin = d;
                first = false;
              }
            }
          }
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          if (z0 >= jz0 && z1 <= jz1) {
            process_.push_back(index);
          } else {
            f |= af::DontIntegrate;
            ignored_.push_back(index);
          }
        } else {
          ignored_.push_back(index);
        }
      }
    }

    af::reflection_table data_;
    vec2<int> array_range_;
    double block_size_;
    std::size_t num_tasks_;
    bool finished_;
    af::shared< tiny<int,2> > jobs_;
    af::shared<std::size_t> ignored_;
    af::shared<std::size_t> process_;
  };


  /**
   * A class to managing spliting and mergin data
   */
  class IntegrationManager3DMultiExecutor {
  public:

    IntegrationManager3DMultiExecutor(
        af::reflection_table reflections,
        vec2<int> array_range,
        double block_size)
          : data_(reflections) {

      // Compute the list of jobs
      jobs_ = detail::compute_jobs(array_range, block_size);

      // Set all the finished flags to false
      finished_.assign(jobs_.size(), false);

      // Generate indices of reflections to be integrated, used as reference
      // spots or passed just as data for each data block. If the reflection is
      // not to be integrated, it is added to each block which it overlaps. If
      // the reflection is a reference spot, it is added to each block in which
      // it is fully recorded. If the spot is to be integrated, it is added to
      // the block in which it is closest to the centre. If the reflection is
      // larger than block_size / 2, then it is not fully recorded in any block
      // and is unprocessed.
      typedef std::pair<std::size_t, bool> job_type;
      typedef std::vector<job_type> job_list_type;

      // Get some reflection data
      af::const_ref<int6> bbox = data_["bbox"];
      af::ref<std::size_t> flags = data_["flags"];

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t jmin = 0;
          double dmin = 0;
          bool first = true;
          for (std::size_t j = 0; j < jobs_.size(); ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              if (f & af::ReferenceSpot) {
                indices[j].push_back(job_type(index, false));
              }
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (first || d < dmin) {
                jmin = j;
                dmin = d;
                first = false;
              }
            }
          }
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          if (z0 >= jz0 && z1 <= jz1) {
            if (f & af::ReferenceSpot) {
              DIALS_ASSERT(indices[jmin].size() > 0);
              DIALS_ASSERT(indices[jmin].back().first == index);
              indices[jmin].back().second = true;
            } else {
              indices[jmin].push_back(job_type(index, true));
            }
          } else {
            f |= af::DontIntegrate;
            ignored_.push_back(index);
          }
        } else {
          ignored_.push_back(index);
        }
      }

      // Compute number of reflections in each task
      std::vector<std::size_t> num(jobs_.size(), 0);
      for (std::size_t i = 0; i < num.size(); ++i) {
        num[i] = indices[i].size();
      }

      // Compute offsests
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));

      // Compute indices
      indices_.resize(offset_.back());
      mask_.resize(offset_.back());
      std::size_t k = 0;
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const job_list_type& ind = indices[i];
        for (std::size_t j = 0; j < ind.size(); ++j, ++k) {
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = ind[j].first;
          mask_[k] = ind[j].second;
          DIALS_ASSERT(indices_[k] < bbox.size());
        }
      }
      DIALS_ASSERT(k == indices_.size());
    }

    /**
     * @returns The result data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return data_;
    }

    /**
     * @returns Is the process finished
     */
    bool finished() const {
      return finished_.all_eq(true);
    }

    /**
     * @returns The number of tasks
     */
    std::size_t size() const {
      return finished_.size();
    }

    /**
     * @returns The block indices
     */
    vec2<int> job(std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The list of reflections to not process.
     */
    af::shared<std::size_t> ignored() const {
      return af::shared<std::size_t>(&ignored_[0], &ignored_[0]+ignored_.size());
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());

      // Extract the reflection table
      af::reflection_table result = select_rows_index(data_, ind);

      // Extract the flags and set those reflections that are not to be
      // processed.
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (msk[i] == false) {
          flags[i] |= af::DontIntegrate;
        }
      }

      // Return the reflections
      return result;
    }

    /**
     * Accumulate the results.
     */
    void accumulate(std::size_t index, af::reflection_table result) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      DIALS_ASSERT(finished_[index] == false);
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());
      DIALS_ASSERT(ind.size() == result.size());

      // Set the result
      set_selected_rows_index_mask(data_, ind, msk, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Get the indices for each job
     */
    af::const_ref<std::size_t> indices(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t> (&indices_[off], num);
    }

    /**
     * Get the mask for each job
     */
    af::const_ref<bool> mask(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= mask_.size());
      return af::const_ref<bool> (&mask_[off], num);
    }

    af::shared< tiny<int,2> > jobs_;
    af::shared<bool> finished_;
    af::reflection_table data_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
    af::shared<bool> mask_;
   af::shared<std::size_t> ignored_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
