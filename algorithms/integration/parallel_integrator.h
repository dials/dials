/*
 * parallel_integrator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H

#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/imageset.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/reflection.h>
#include <dials/error.h>
#include <dials/util/thread_pool.h>
#include <dials/algorithms/shoebox/find_overlapping.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/centroid/centroid.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <map>

#include <dials/algorithms/integration/interfaces.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dxtbx::model::Panel;

  using dxtbx::ImageSweep;
  using dxtbx::format::Image;
  using dxtbx::format::ImageTile;

  using dials::model::Shoebox;
  using dials::model::AdjacencyList;

  /**
   * Class to wrap logging
   */
  class Logger {
  public:

    Logger(boost::python::object obj)
      : obj_(obj) {}

    void info(const char *str) const {
      obj_.attr("info")(str);
    }

    void debuf(const char *str) const {
      obj_.attr("debug")(str);
    }

  private:

    boost::python::object obj_;

  };


  /**
   * A class to store the image data buffer
   */
  class Buffer {
  public:

    /**
     * Initialise the the size of the panels
     * @param detector The detector model
     * @param num_images The number of images
     * @param use_dynamic_mask Using a dynamic mask
     * @param external_mask The external mask
     */
    Buffer(const Detector &detector,
           std::size_t num_images,
           bool use_dynamic_mask,
           const Image<bool> &external_mask)
        : use_dynamic_mask_(use_dynamic_mask) {
      std::size_t zsize = num_images;
      DIALS_ASSERT(zsize > 0);
      for (std::size_t i = 0; i < detector.size(); ++i) {
        std::size_t xsize = detector[i].get_image_size()[0];
        std::size_t ysize = detector[i].get_image_size()[1];
        DIALS_ASSERT(xsize > 0);
        DIALS_ASSERT(ysize > 0);

        // Allocate all the data buffers
        data_.push_back(
            af::versa< double, af::c_grid<3> >(
              af::c_grid<3>(zsize, ysize, xsize)));

        // Allocate all the dynamic mask buffers
        if (use_dynamic_mask) {
          dynamic_mask_.push_back(
            af::versa< bool, af::c_grid<3> >(
              af::c_grid<3>(zsize, ysize, xsize)));
        }

        // Allocate the static mask buffer
        static_mask_.push_back(
            af::versa< bool, af::c_grid<2> >(
              af::c_grid<2>(ysize, xsize), true));
      }

      // Set the external mask
      if (!external_mask.empty()) {
        DIALS_ASSERT(external_mask.n_tiles() == static_mask_.size());
        for (std::size_t i = 0; i < external_mask.n_tiles(); ++i) {
          set_external_mask_for_panel(
              external_mask.tile(i).data().const_ref(),
              static_mask_[i].ref());
        }
      }
    }

    /**
     * Copy an image to the buffer
     * @param data The image data
     * @param index The image index
     */
    void copy(const Image<double> &data, std::size_t index) {
      DIALS_ASSERT(use_dynamic_mask_ == false);
      DIALS_ASSERT(data.n_tiles() == data_.size());
      for (std::size_t i = 0; i < data.n_tiles(); ++i) {
        copy(data.tile(i).data().const_ref(), data_[i].ref(), index);
      }
    }

    /**
     * Copy an image to the buffer
     * @param data The image data
     * @param mask The mask data
     * @param index The image index
     */
    void copy(const Image<double> &data, const Image<bool> &mask, std::size_t index) {
      DIALS_ASSERT(use_dynamic_mask_ == true);
      DIALS_ASSERT(data.n_tiles() == mask.n_tiles());
      DIALS_ASSERT(data.n_tiles() == data_.size());
      DIALS_ASSERT(mask.n_tiles() == dynamic_mask_.size());
      for (std::size_t i = 0; i < data.n_tiles(); ++i) {
        copy(data.tile(i).data().const_ref(), data_[i].ref(), index);
        copy(mask.tile(i).data().const_ref(), dynamic_mask_[i].ref(), index);
        apply_static_mask(dynamic_mask_[i].ref(), static_mask_[i].const_ref(), index);
      }
    }

    /**
     * @param The panel number
     * @returns The buffer for the panel
     */
    af::const_ref< double, af::c_grid<3> > data(std::size_t panel) const {
      DIALS_ASSERT(panel < data_.size());
      return data_[panel].const_ref();
    }

    /**
     * @param The panel number
     * @returns The buffer for the panel
     */
    af::const_ref< bool, af::c_grid<3> > dynamic_mask(std::size_t panel) const {
      DIALS_ASSERT(use_dynamic_mask_ == true);
      DIALS_ASSERT(panel < dynamic_mask_.size());
      return dynamic_mask_[panel].const_ref();
    }

    /**
     * @param The panel number
     * @returns The buffer for the panel
     */
    af::const_ref< bool, af::c_grid<2> > static_mask(std::size_t panel) const {
      DIALS_ASSERT(use_dynamic_mask_ == false);
      DIALS_ASSERT(panel < static_mask_.size());
      return static_mask_[panel].const_ref();
    }

    /**
     * @returns do we have a dynamic mask?
     */
    bool has_dynamic_mask() const {
      return use_dynamic_mask_;
    }

  protected:

    /**
     * Copy the data from 1 panel
     * @param src The source
     * @param dst The destination
     * @param index The image index
     */
    template <typename T>
    void copy(af::const_ref< T, af::c_grid<2> > src,
              af::ref < T, af::c_grid<3> > dst,
              std::size_t index) {
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];
      DIALS_ASSERT(index < dst.accessor()[0]);
      DIALS_ASSERT(src.accessor()[0] == dst.accessor()[1]);
      DIALS_ASSERT(src.accessor()[1] == dst.accessor()[2]);
      for (std::size_t j = 0; j < ysize * xsize; ++j) {
        dst[index * (xsize*ysize) + j] = src[j];
      }
    }

    /**
     * Apply the static mask to the dynamic mask
     * @param src The source
     * @param dst The destination
     * @param index The image index
     */
    void apply_static_mask(af::ref< bool, af::c_grid<3> > dynamic_mask,
                           af::const_ref< bool, af::c_grid<2> > static_mask,
                           std::size_t index) {
      std::size_t ysize = dynamic_mask.accessor()[1];
      std::size_t xsize = dynamic_mask.accessor()[2];
      DIALS_ASSERT(index < dynamic_mask.accessor()[0]);
      DIALS_ASSERT(dynamic_mask.accessor()[1] == static_mask.accessor()[0]);
      DIALS_ASSERT(dynamic_mask.accessor()[2] == static_mask.accessor()[1]);
      for (std::size_t j = 0; j < ysize * xsize; ++j) {
        bool d = dynamic_mask[index * (xsize*ysize) + j];
        bool s = static_mask[j];
        dynamic_mask[index * (xsize*ysize) + j] = d && s;
      }
    }

    /**
     * Set the static mask for a panel
     * @param panel The panel
     * @param panel_static_mask The static mask
     */
    void set_external_mask_for_panel(
            const af::const_ref< bool, af::c_grid<2> > &external_mask,
            af::ref< bool, af::c_grid<2> > panel_static_mask) {
      DIALS_ASSERT(external_mask.accessor().all_eq(panel_static_mask.accessor()));
      for (std::size_t j = 0; j < external_mask.size(); ++j) {
        panel_static_mask[j] = external_mask[j] && panel_static_mask[j];
      }
    }

    std::vector< af::versa< double, af::c_grid<3> > > data_;
    std::vector< af::versa< bool, af::c_grid<3> > > dynamic_mask_;
    std::vector< af::versa< bool, af::c_grid<2> > > static_mask_;
    bool use_dynamic_mask_;

  };


  /**
   * A class to integrate a single reflection
   */
  class ReflectionIntegrator {
  public:

    /**
     * Initialise the integrator
     * @param compute_mask The mask calculation function
     * @param compute_background The background calculation function
     * @param compute_intensity The intensity calculation function
     * @param buffer The image buffer array
     * @param zstart The first image index
     * @param underload The underload value
     * @param overload The overload value
     */
    ReflectionIntegrator(
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity,
          const Buffer &buffer,
          int zstart,
          double underload,
          double overload,
          bool debug)
      : compute_mask_(compute_mask),
        compute_background_(compute_background),
        compute_intensity_(compute_intensity),
        buffer_(buffer),
        zstart_(zstart),
        underload_(underload),
        overload_(overload),
        debug_(debug) {}

    /**
     * Integrate a reflection using the following procedure:
     *
     * 1. Extract the buffered image data into a shoebox
     * 2. Compute the mask for the reflection
     * 3. Compute the mask for adjacent reflections
     * 4. Compute the reflection background
     * 5. Compute the reflection centroid
     * 6. Compute the summed intensity
     * 7. Compute the profile fitted intensity
     * 8. Delete the shoebox unless debuf has been set
     *
     * @param reflection The reflection object
     * @param adjacent_reflections The list of adjacent reflections
     */
    void operator()(af::Reflection &reflection,
                    std::vector<af::Reflection> &adjacent_reflections) const {

      // Get the panel number
      std::size_t panel = reflection.get<std::size_t>("panel");

      // Extract the shoebox data
      if (buffer_.has_dynamic_mask()) {
        extract_shoebox(
            buffer_.data(panel),
            buffer_.dynamic_mask(panel),
            reflection,
            zstart_,
            underload_,
            overload_);
      } else {
        extract_shoebox(
            buffer_.data(panel),
            buffer_.static_mask(panel),
            reflection,
            zstart_,
            underload_,
            overload_);
      }

      // Compute the mask
      compute_mask_(reflection);

      // Set all the bounding boxes of adjacent reflections
      // And compute the mask for these reflections too.
      for (std::size_t i = 0; i < adjacent_reflections.size(); ++i) {
        adjacent_reflections[i]["bbox"] = reflection.get<int6>("bbox");
        adjacent_reflections[i]["shoebox"] = reflection.get< Shoebox<> >("shoebox");
        compute_mask_(adjacent_reflections[i], true);
      }

      // Compute the background
      try {
        compute_background_(reflection);
      } catch (dials::error) {
        delete_shoebox(reflection, adjacent_reflections);
        return;
      }

      // Compute the centroid
      compute_centroid(reflection);

      // Compute the summed intensity
      compute_summed_intensity(reflection);

      // Compute the profile fitted intensity
      try {
        compute_intensity_(reflection, adjacent_reflections);
      } catch (dials::error) {
        // pass
      }

      // Erase the shoebox
      delete_shoebox(reflection, adjacent_reflections);
    }


  protected:

    /**
     * Delete the shoebox
     * @param reflection The reflection
     * @param adjacent_reflections The adjancent reflections
     */
    void delete_shoebox(
          af::Reflection &reflection,
          std::vector<af::Reflection> &adjacent_reflections) const {

      // Erase the shoebox from the reflection
      if (!debug_) {
        reflection.erase("shoebox");
        for (std::size_t i = 0; i < adjacent_reflections.size(); ++i) {
          adjacent_reflections[i].erase("shoebox");
        }
      }
    }

    /**
     * Extract the shoebox data from the buffer
     */
    void extract_shoebox(
          const af::const_ref< double, af::c_grid<3> > &data_buffer,
          const af::const_ref< bool, af::c_grid<2> > &mask_buffer,
          af::Reflection &reflection,
          int zstart,
          double underload,
          double overload) const {
      DIALS_ASSERT(data_buffer.accessor()[1] == mask_buffer.accessor()[0]);
      DIALS_ASSERT(data_buffer.accessor()[2] == mask_buffer.accessor()[1]);
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref< float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref< int,   af::c_grid<3> > mask = shoebox.mask.ref();
      int x0 = bbox[0];
      int x1 = bbox[1];
      int y0 = bbox[2];
      int y1 = bbox[3];
      int z0 = bbox[4];
      int z1 = bbox[5];
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(z1 > z0);
      std::size_t zsize = z1 - z0;
      std::size_t ysize = y1 - y0;
      std::size_t xsize = x1 - x0;
      DIALS_ASSERT(zsize == data.accessor()[0]);
      DIALS_ASSERT(ysize == data.accessor()[1]);
      DIALS_ASSERT(xsize == data.accessor()[2]);
      DIALS_ASSERT(shoebox.is_consistent());
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int kk = z0 + k - zstart;
            int jj = y0 + j;
            int ii = x0 + i;
            if (kk >= 0 &&
                jj >= 0 &&
                ii >= 0 &&
                kk < data_buffer.accessor()[0] &&
                jj < data_buffer.accessor()[1] &&
                ii < data_buffer.accessor()[2]) {
              double d = data_buffer(kk, jj, ii);
              int m = mask_buffer(jj, ii) && (d > underload && d < overload)
                ? Valid
                : 0;
              data(k,j,i) = d;
              mask(k,j,i) = m;
            } else {
              data(k,j,i) = 0;
              mask(k,j,i) = 0;
            }
          }
        }
      }
      reflection["shoebox"] = shoebox;
    }

    /**
     * Extract the shoebox data and mask from the buffer
     */
    void extract_shoebox(
          const af::const_ref< double, af::c_grid<3> > &data_buffer,
          const af::const_ref< bool, af::c_grid<3> > &mask_buffer,
          af::Reflection &reflection,
          int zstart,
          double underload,
          double overload) const {
      DIALS_ASSERT(data_buffer.accessor().all_eq(mask_buffer.accessor()));
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref< float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref< int,   af::c_grid<3> > mask = shoebox.mask.ref();
      int x0 = bbox[0];
      int x1 = bbox[1];
      int y0 = bbox[2];
      int y1 = bbox[3];
      int z0 = bbox[4];
      int z1 = bbox[5];
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(z1 > z0);
      std::size_t zsize = z1 - z0;
      std::size_t ysize = y1 - y0;
      std::size_t xsize = x1 - x0;
      DIALS_ASSERT(zsize == data.accessor()[0]);
      DIALS_ASSERT(ysize == data.accessor()[1]);
      DIALS_ASSERT(xsize == data.accessor()[2]);
      DIALS_ASSERT(shoebox.is_consistent());
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int kk = z0 + k - zstart;
            int jj = y0 + j;
            int ii = x0 + i;
            if (kk >= 0 &&
                jj >= 0 &&
                ii >= 0 &&
                kk < data_buffer.accessor()[0] &&
                jj < data_buffer.accessor()[1] &&
                ii < data_buffer.accessor()[2]) {
              double d = data_buffer(kk, jj, ii);
              int m = mask_buffer(kk, jj, ii) && (d > underload && d < overload)
                ? Valid
                : 0;
              data(k,j,i) = d;
              mask(k,j,i) = m;
            } else {
              data(k,j,i) = 0;
              mask(k,j,i) = 0;
            }
          }
        }
      }
      reflection["shoebox"] = shoebox;
    }

    /**
     * Compute the centroid
     */
    void compute_centroid(af::Reflection &reflection) const {

      using dials::model::Centroid;

      // Get the shoebox and compute centroid
      Shoebox<> shoebox = reflection.get< Shoebox<> >("shoebox");
      Centroid centroid = shoebox.centroid_foreground_minus_background();

      // Set the centroid values
      reflection["xyzobs.px.value"] = centroid.px.position;
      reflection["xyzobs.px.variance"] = centroid.px.variance;
    }

    /**
     * Compute the summed intensity
     */
    void compute_summed_intensity(af::Reflection &reflection) const {

      using dials::model::Intensity;

      // Get flags and reset
      std::size_t flags = reflection.get<std::size_t>("flags");
      flags &= ~af::IntegratedSum;
      flags &= ~af::FailedDuringSummation;

      // Get the shoebox and compute the summed intensity
      Shoebox<> shoebox = reflection.get< Shoebox<> >("shoebox");
      Intensity intensity = shoebox.summed_intensity();

      // Set the intensities
      reflection["intensity.sum.value"] = intensity.observed.value;
      reflection["intensity.sum.variance"] = intensity.observed.variance;
      reflection["background.sum.value"] = intensity.background.value;
      reflection["background.sum.variance"] = intensity.background.variance;

      // Set the appropriate flag
      if (intensity.observed.success) {
        flags |= af::IntegratedSum;
      } else {
        flags |= af::FailedDuringSummation;
      }
      reflection["flags"] = flags;
    }

    const MaskCalculatorIface &compute_mask_;
    const BackgroundCalculatorIface &compute_background_;
    const IntensityCalculatorIface &compute_intensity_;
    const Buffer &buffer_;
    int zstart_;
    double underload_;
    double overload_;
    bool debug_;
  };


  /**
   * a class to sort the indices of all reflections that are fully recorded
   * after a particular image.
   */
  class Lookup {
  public:

    /**
     * @param bbox the bounding boxs
     * @param zstart the first frame number
     * @param n the number of frames
     */
    Lookup(af::const_ref<int6> bbox, int zstart, std::size_t n)
      : indices_(bbox.size()) {

      // fill the index array
      for (std::size_t i = 0; i < indices_.size(); ++i) {
        indices_[i] = i;
      }

      // sort the indices by the final frame number
      std::sort(indices_.begin(), indices_.end(), sort_by_frame(bbox));
      DIALS_ASSERT(bbox[indices_.front()][5] - zstart >= 1);
      DIALS_ASSERT(bbox[indices_.back()][5]  - zstart <= n);

      // create an offset array that records the positions in the index array
      // where the frame increments such that
      // offset[i], offset[i+1] gives the range of indices for a frame.
      std::size_t i = 0;
      offset_.push_back(0);
      for (std::size_t j = 0; j < n; ++j) {
        while (i < indices_.size() && bbox[indices_[i]][5] - zstart <= j+1) i++;
        offset_.push_back(i);
      }
      DIALS_ASSERT(offset_.size() == n + 1);
      DIALS_ASSERT(offset_.back() == indices_.size());
    }

    /**
     * @param z the frame number
     * @returns the indices for a given frame
     */
    af::const_ref<std::size_t> indices(std::size_t z) const {
      DIALS_ASSERT(z < offset_.size()-1);
      DIALS_ASSERT(offset_[z+1] >= offset_[z]);
      std::size_t i = offset_[z];
      std::size_t n = offset_[z+1] - i;
      return af::const_ref<std::size_t>(&indices_[i], n);
    }

  private:

    /**
     * helper function to sort by final bbox frame
     */
    struct sort_by_frame {
      af::const_ref<int6> bbox_;
      sort_by_frame(af::const_ref<int6> bbox) : bbox_(bbox) {}
      bool operator()(std::size_t a, std::size_t b) const {
        return bbox_[a][5] < bbox_[b][5];
      }
    };

    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
  };



  /**
   * A class to perform parallel integration
   */
  class ParallelIntegrator {
  public:

    /**
     * Do the integration
     * @param reflections The reflection table
     * @param imageset The imageset
     * @param compute_mask The mask calulcation function
     * @param compute_background The background calculation function
     * @param compute_intensity The intensity calculation function
     * @param logger The logger class
     * @param nthreads The number of parallel threads
     * @param use_dynamic_mask Use the dynamic mask if present
     * @param debug Add debug output
     */
    ParallelIntegrator(
          af::reflection_table reflections,
          ImageSweep imageset,
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity,
          const Logger &logger,
          std::size_t nthreads,
          bool use_dynamic_mask,
          bool debug) {

      using dials::algorithms::shoebox::find_overlapping_multi_panel;

      // Check the input
      DIALS_ASSERT(nthreads > 0);

      // Check the models
      DIALS_ASSERT(imageset.get_detector() != NULL);
      DIALS_ASSERT(imageset.get_scan() != NULL);
      Detector detector = *imageset.get_detector();
      Scan scan = *imageset.get_scan();

      // Get the size of the data buffer needed
      std::size_t zsize = imageset.size();
      DIALS_ASSERT(zsize > 0);

      // Get the starting frame and the underload/overload values
      int zstart = scan.get_array_range()[0];
      double underload = detector[0].get_trusted_range()[0];
      double overload  = detector[0].get_trusted_range()[1];
      DIALS_ASSERT(underload < overload);
      for (std::size_t i = 1; i < detector.size(); ++i) {
        DIALS_ASSERT(underload == detector[i].get_trusted_range()[0]);
        DIALS_ASSERT(overload  == detector[i].get_trusted_range()[1]);
      }

      // Get the reflection flags and bbox
      af::const_ref<std::size_t> panel = reflections.get<std::size_t>("panel").const_ref();
      af::const_ref<int6> bbox = reflections.get<int6>("bbox").const_ref();
      af::ref<std::size_t> flags = reflections.get<std::size_t>("flags").ref();

      // Reset the flags
      reset_flags(flags);

      // Find the overlapping reflections
      AdjacencyList overlaps = find_overlapping_multi_panel(bbox, panel);

      // Allocate the array for the image data
      Buffer buffer(
          detector,
          zsize,
          use_dynamic_mask && imageset.has_dynamic_mask(),
          imageset.get_static_mask());

      // If we have shoeboxes then delete
      if (reflections.contains("shoebox")) {
        reflections.erase("shoebox");
      }

      // Transform reflection data from column major to row major. The
      // reason for doing this is because we want to process each reflection
      // in parallel. Therefore passing a reflection object is better than
      // passing the whole table. Additionally, because the reflection table
      // uses a std::map, it may not be thread safe, so we want to avoid
      // accessing this across multiple threads.
      af::shared<af::Reflection> reflection_array =
        reflection_table_to_array(reflections);

      // The lookup class gives the indices of reflections whose bounding boxes
      // are complete at a given image. This is used to submit reflections for
      // integration after each image is processed.
      Lookup lookup(bbox, zstart, zsize);

      // Create the reflection integrator. This class is called for each
      // reflection to integrate the data
      ReflectionIntegrator integrator(
          compute_mask,
          compute_background,
          compute_intensity,
          buffer,
          zstart,
          underload,
          overload,
          debug);

      // Do the integration
      process(
          lookup,
          integrator,
          buffer,
          reflection_array.ref(),
          overlaps,
          imageset,
          bbox,
          flags,
          nthreads,
          logger);

      // Transform the row major reflection array to the reflection table
      reflections_ = reflection_table_from_array(reflection_array.const_ref());
    }

    /**
     * @returns The integrated reflections
     */
    af::reflection_table reflections() const {
      return reflections_;
    }

    /**
     * Static method to get the memory in bytes needed
     * @param imageset the imageset class
     * @param use_dynamic_mask Are we using dynamic mask
     */
    static
    std::size_t compute_required_memory(
        ImageSweep imageset,
        bool use_dynamic_mask) {
      DIALS_ASSERT(imageset.get_detector() != NULL);
      DIALS_ASSERT(imageset.get_scan() != NULL);
      Detector detector = *imageset.get_detector();
      Scan scan = *imageset.get_scan();
      std::size_t nelements = 0;
      for (std::size_t i = 0; i < detector.size(); ++i) {
        std::size_t xsize = detector[i].get_image_size()[0];
        std::size_t ysize = detector[i].get_image_size()[1];
        nelements += xsize * ysize;
      }
      nelements *= scan.get_num_images();
      std::size_t nbytes = nelements * sizeof(double);
      if (use_dynamic_mask) {
        nbytes += nelements * sizeof(bool);
      }
      return nbytes;
    }

    /**
     * Static method to get the memory in bytes needed
     * @param imageset the imageset class
     * @param use_dynamic_mask Are we using dynamic mask
     * @param max_memory_usage The maximum memory usage
     */
    static
    std::size_t compute_max_block_size(
        ImageSweep imageset,
        bool use_dynamic_mask,
        std::size_t max_memory_usage) {
      DIALS_ASSERT(max_memory_usage > 0);
      DIALS_ASSERT(imageset.get_detector() != NULL);
      Detector detector = *imageset.get_detector();
      std::size_t nelements = 0;
      for (std::size_t i = 0; i < detector.size(); ++i) {
        std::size_t xsize = detector[i].get_image_size()[0];
        std::size_t ysize = detector[i].get_image_size()[1];
        nelements += xsize * ysize;
      }
      std::size_t nbytes = nelements * sizeof(double);
      if (use_dynamic_mask) {
        nbytes += nelements * sizeof(bool);
      }
      DIALS_ASSERT(nbytes > 0);
      DIALS_ASSERT(max_memory_usage > nbytes);
      return (std::size_t)std::floor((float)max_memory_usage / (float)nbytes);
    }

  protected:

    /**
     * Reset the reflection flags
     */
    void reset_flags(af::ref<std::size_t> flags) const {
      for (std::size_t i = 0; i < flags.size(); ++i) {
        flags[i] &= ~af::IntegratedSum;
        flags[i] &= ~af::IntegratedPrf;
      }
    }

    /**
     * Do the processing by the following procedure.
     *
     * 1. For each reflection, construct a list of adjacent reflection
     * 2. Loop through all the images in the imageset
     * 3. Copy the images to a buffer
     * 4. Loop through all reflections complete on the image
     * 5. For each complete reflection post a reflection integration job to the
     *    thread pool.
     */
    void process(
        const Lookup &lookup,
        const ReflectionIntegrator &integrator,
        Buffer &buffer,
        af::ref< af::Reflection > reflections,
        const AdjacencyList &overlaps,
        ImageSweep imageset,
        af::const_ref<int6> bbox,
        af::const_ref<std::size_t> flags,
        std::size_t nthreads,
        const Logger &logger) const {

      using dials::util::ThreadPool;

      // Construct a list of adjacent reflections for each given reflection
      std::vector< std::vector<af::Reflection> > adjacent_reflections(reflections.size());
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if ((flags[i] & af::DontIntegrate) == 0) {
          adjacent_reflections[i].reserve(overlaps.vertex_num_edges(i));
          AdjacencyList::edge_iterator_range edges = overlaps.edges(i);
          for (AdjacencyList::edge_iterator it = edges.first; it != edges.second; ++it) {
            DIALS_ASSERT(it->first == i);
            DIALS_ASSERT(it->second < reflections.size());
            adjacent_reflections[i].push_back(reflections[it->second]);
          }
        }
      }

      // Create the thread pool
      ThreadPool pool(nthreads);

      // Get the size of the array
      int zstart = imageset.get_scan()->get_array_range()[0];
      std::size_t zsize = imageset.size();

      // Loop through all the images
      for (std::size_t i = 0; i < zsize; ++i) {

        // Copy the image to the buffer
        if (buffer.has_dynamic_mask()) {
          buffer.copy(imageset.get_corrected_data(i), imageset.get_dynamic_mask(i), i);
        } else {
          buffer.copy(imageset.get_corrected_data(i), i);
        }

        // Get the reflections recorded at this point
        af::const_ref<std::size_t> indices = lookup.indices(i);

        // Iterate through the reflection indices
        std::size_t count = 0;
        for (std::size_t j = 0; j < indices.size(); ++j) {

          // Get the reflection index
          std::size_t k = indices[j];

          // Check that the reflection bounding box is within the range of
          // images that have been read
          DIALS_ASSERT(bbox[k][5] <= zstart + i + 1);

          // Ignore if we're not integrating this reflection
          if (flags[k] & af::DontIntegrate) {
            continue;
          } else {
            count++;
          }

          // Post the integration job
          pool.post(
            boost::bind(
              &ReflectionIntegrator::operator(),
              integrator,
              boost::ref(reflections[k]),
              boost::ref(adjacent_reflections[k])));
        }

        // Print some output
        std::ostringstream ss;
        ss << "Integrating "
           << count
           << " reflections on image "
           << zstart + i;
        logger.info(ss.str().c_str());
      }

      // Wait for all the integration jobs to complete
      pool.wait();
    }

    af::reflection_table reflections_;
  };


  /**
   * A class to manage jobs
   */
  class SimpleJobList {
  public:


    /**
     * Compute the jobs
     * @param range The range of frames
     * @param block_size The size of the blocks
     */
    SimpleJobList(tiny<int,2> range, int block_size) {
      construct_job_list(range, block_size);
      construct_frame_to_job_lookup();
    }

    /**
     * Set the jobs
     * @params jobs The list of jobs
     */
    SimpleJobList(const af::const_ref< tiny<int,2> > &jobs) {
      DIALS_ASSERT(jobs.size() > 0);
      DIALS_ASSERT(jobs[0][1] > jobs[0][0]);
      jobs_.push_back(jobs[0]);
      for (std::size_t i = 1; i < jobs.size(); ++i) {
        DIALS_ASSERT(jobs[i][1] > jobs[i][0]);
        DIALS_ASSERT(jobs[i][0] > jobs[i-1][0]);
        DIALS_ASSERT(jobs[i][1] > jobs[i-1][1]);
        DIALS_ASSERT(jobs[i][0] <= jobs[i-1][1]);
        jobs_.push_back(jobs[i]);
      }
      construct_frame_to_job_lookup();
    }

    /**
     * @returns The requested job
     */
    tiny<int,2> operator[](std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The number of jobs
     */
    std::size_t size() const {
      return jobs_.size();
    }

    /**
     * Get the index of the job closest to this frame
     * @param frame The frame number
     * @returns The job index
     */
    std::size_t job_index(int frame) const {
      std::size_t index = frame - jobs_.front()[0];
      DIALS_ASSERT(index >= 0);
      DIALS_ASSERT(index < frame_to_job_lookup_.size());
      return frame_to_job_lookup_[index];
    }

  private:

    /**
     * Construct the job list
     * @param range The range of frames
     * @param block_size The block size
     */
    void construct_job_list(tiny<int,2> range, int block_size) {

      // Check some input
      int frame0 = range[0];
      int frame1 = range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      DIALS_ASSERT(nframes > 0);

      // Block size is clamped to number of frames
      if (block_size > nframes) {
        block_size = nframes;
      }
      DIALS_ASSERT(block_size > 0);

      // If the block size is equal to 1, then add all frames as jobs
      // otherwise compute the jobs
      if (block_size == 1) {
        for (int f = frame0; f < frame1; ++f) {
          jobs_.push_back(tiny<int,2>(f, f+1));
        }
      } else {

        // Compute the half block size such that images are divided between
        // blocks more evenly spaced
        int nblocks = (int)std::ceil(2.0 * nframes / (double)block_size);
        DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
        int half_block_size = (int)std::ceil((double)nframes / (double)nblocks);

        // Construct a list of job indices
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

        // Add all the jobs to the list
        DIALS_ASSERT(indices.front() == frame0);
        DIALS_ASSERT(indices.back() == frame1);
        DIALS_ASSERT(indices.size() > 2);
        for (std::size_t i = 0; i < indices.size() - 2; ++i) {
          int i1 = indices[i];
          int i2 = indices[i+2];
          DIALS_ASSERT(i2 > i1);
          jobs_.push_back(tiny<int,2>(i1, i2));
        }
        DIALS_ASSERT(jobs_.size() > 0);
      }
    }

    /**
     * Construct the frame to job lookup table
     */
    void construct_frame_to_job_lookup() {

      // Check all the jobs overlap and are in order
      for (std::size_t i = 0; i < jobs_.size()-1; ++i) {
        DIALS_ASSERT(jobs_[i][0] < jobs_[i][1]);
        DIALS_ASSERT(jobs_[i+1][0] < jobs_[i+1][1]);
        DIALS_ASSERT(jobs_[i][0] < jobs_[i+1][0]);
        DIALS_ASSERT(jobs_[i][1] >= jobs_[i+1][0]);
        DIALS_ASSERT(jobs_[i][1] < jobs_[i+1][1]);
      }

      // set the first and last frames
      int first_frame = jobs_.front()[0];
      int last_frame = jobs_.back()[1];
      DIALS_ASSERT(first_frame < last_frame);

      // Loop through all the frames and find the job to which the frame is
      // closest to the centre. Add the frame to that job
      std::size_t closest_index = 0;
      for (int frame = first_frame; frame < last_frame; ++frame) {
        int z0 = jobs_[closest_index][0];
        int z1 = jobs_[closest_index][1];
        double zc = (z0 + z1) / 2.0;
        double closest_distance = std::abs(zc - (frame + 0.5));
        for (std::size_t i = closest_index+1; i < jobs_.size(); ++i) {
          int zz0 = jobs_[i][0];
          int zz1 = jobs_[i][1];
          double zzc = (zz0 + zz1) / 2.0;
          double distance = std::abs(zzc - (frame + 0.5));
          if (distance < closest_distance) {
            closest_distance = distance;
            closest_index = i;
          } else {
            break;
          }
        }
        frame_to_job_lookup_.push_back(closest_index);
      }

    }

    std::vector< tiny<int,2> > jobs_;
    std::vector< std::size_t > frame_to_job_lookup_;
  };




  /**
   * A lookup class to split reflections along job boundaries and to gives
   * indices of reflections in each job
   */
  class SimpleReflectionLookup {
  public:

    /**
     * Construct the lookup
     * @param jobs The job list
     * @param data The reflections
     */
    SimpleReflectionLookup(
          const SimpleJobList &jobs,
          af::reflection_table data)
        : jobs_(jobs),
          first_frame_(0),
          last_frame_(0) {

      // Split the reflections along job boundaries
      data_ = split_reflections(data);

      // Construct the lookup table
      construct_job_to_reflection_lookup();
    }

    /**
     * @returns The reflection data
     */
    af::reflection_table data() {
      return data_;
    }

    /**
     * Get the indices of reflections in this job
     * @param index The job index
     * @returns The list of reflection indices
     */
    af::const_ref<std::size_t> indices(std::size_t index) const {
      DIALS_ASSERT(index < job_to_reflection_lookup_.size());
      return af::const_ref<std::size_t>(
          &job_to_reflection_lookup_[index][0],
          job_to_reflection_lookup_[index].size());
    }

    /**
     * Get the index of the job closest to this frame
     * @param frame The frame number
     * @returns The job index
     */
    std::size_t job_index(int frame) const {
      return jobs_.job_index(frame);
    }

    /**
     * Get the frame range for the job
     * @param index The job index
     * @returns The frames in the job
     */
    tiny<int,2> job_range(std::size_t index) const {
      return jobs_[index];
    }

  protected:

    /**
     * Create the lookup of reflections in each job
     */
    void construct_job_to_reflection_lookup() {
      DIALS_ASSERT(data_.is_consistent());
      DIALS_ASSERT(data_.size() > 0);
      DIALS_ASSERT(data_.contains("bbox"));
      DIALS_ASSERT(jobs_.size() > 0);

      // Get some arrays
      af::const_ref<int6> bbox = data_["bbox"];

      // For each reflection, get the frame range and then get the centre of
      // that range. Now find the job closest to that z point. Check that the
      // splitting has worked correctly and that the reflection is within the
      // job frame range. Then add the reflection to the list for the block.
      job_to_reflection_lookup_.resize(jobs_.size());
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        int z0 = bbox[i][4];
        int z1 = bbox[i][5];
        int zc = (int)std::floor((z0 + z1) / 2.0);
        int index = job_index(zc);
        tiny<int,2> job = job_range(index);
        DIALS_ASSERT(z0 >= job[0]);
        DIALS_ASSERT(z1 <= job[1]);
        DIALS_ASSERT(index < job_to_reflection_lookup_.size());
        job_to_reflection_lookup_[index].push_back(i);
      }
    }

    /**
     * Split the reflections overlapping job boundaries
     * @param data The input reflection table
     * @returns The split reflection table
     */
    af::reflection_table split_reflections(af::reflection_table data) const {

      // Check input
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.contains("bbox"));
      DIALS_ASSERT(jobs_.size() > 0);

      // Get some arrays
      af::const_ref<int6> bbox = data["bbox"];

      // Split the reflection
      af::shared<std::size_t> indices;
      af::shared<int6> bbox_new;
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        int z0 = bbox[i][4];
        int z1 = bbox[i][5];
        DIALS_ASSERT(z0 < z1);
        std::vector< tiny<int,2> > splits;
        split_at_boundaries(z0, z1, std::back_inserter(splits));
        DIALS_ASSERT(splits.size() > 0);
        for (std::size_t j = 0; j < splits.size(); ++j) {
          int6 b = bbox[i];
          b[4] = splits[j][0];
          b[5] = splits[j][1];
          indices.push_back(i);
          bbox_new.push_back(b);
        }
      }

      // Resize the reflection table
      DIALS_ASSERT(bbox_new.size() == indices.size());
      DIALS_ASSERT(bbox_new.size() >= bbox.size());
      data.resize(bbox_new.size());

      // Reorder the reflections
      af::boost_python::flex_table_suite::reorder(data, indices.const_ref());

      // Set the new bounding boxes
      af::boost_python::flex_table_suite::setitem_column(
          data, "bbox", bbox_new.const_ref());
      af::boost_python::flex_table_suite::setitem_column(
          data, "partial_id", indices.const_ref());

      // Return the data
      return data;
    }

    /**
     * Split a reflection at job boundaries so that most of the reflection is
     * recorded within a single job. Works recursively.
     * @param z0 The first frame in the shoebox
     * @param z1 The last frame in the shoebox
     * @param out The output iterator
     */
    template <typename OutputIterator>
    void split_at_boundaries(int z0, int z1, OutputIterator out) const {
      DIALS_ASSERT(z0 < z1);

      // Compute the job closest to the centre
      int zc = (int)std::floor((z0 + z1) / 2.0);
      int index = job_index(zc);
      tiny<int,2> job = job_range(index);

      // Get the min and max frame range
      int zmin = std::max(z0, job[0]);
      int zmax = std::min(z1, job[1]);
      DIALS_ASSERT(zmin < zmax);

      // If the reflection extends below the job range
      // then split at the lower boundary
      if (z0 < zmin) {
        split_at_boundaries(z0, zmin, out);
      }

      // Append the new frame range
      *out++ = tiny<int,2>(zmin, zmax);

      // If the reflection extends above the job range
      // then split at the upper boundary
      if (z1 > zmax) {
        split_at_boundaries(zmax, z1, out);
      }
    }

    SimpleJobList jobs_;
    int first_frame_;
    int last_frame_;
    af::reflection_table data_;
    std::vector< std::vector<std::size_t> >job_to_reflection_lookup_;
  };


  /**
   * A class to manage the reflections
   */
  class SimpleReflectionManager {
  public:

    /**
     * Construct from the job list and reflection data
     * @param jobs The job list
     * @param data The reflection data
     */
    SimpleReflectionManager(
          const SimpleJobList &jobs,
          af::reflection_table data)
      : lookup_(jobs, data),
        finished_(jobs.size()) {

    }

    /**
     * @returns The result data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return lookup_.data();
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
     * @returns The job
     */
    tiny<int,2> job(std::size_t index) const {
      return lookup_.job_range(index);
    }

    /**
     * @returns The number of reflections in a job
     */
    std::size_t num_reflections(std::size_t index) const {
      DIALS_ASSERT(index < finished_.size());
      return lookup_.indices(index).size();
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {
      using namespace af::boost_python::flex_table_suite;
      DIALS_ASSERT(index < finished_.size());

      // Get the job range
      tiny<int,2> frame = lookup_.job_range(index);

      // Select reflections to process in block
      af::reflection_table data = select_rows_index(
          lookup_.data(),
          lookup_.indices(index));

      // Select other reflections from adjacent blocks. These reflections will
      // not be processed but will be used in finding adjacent reflections.
      // First select reflections from the block before.
      if (index > 0) {
        af::reflection_table temp = select_rows_index(
            lookup_.data(),
            lookup_.indices(index-1));
        af::ref<int6> bbox = temp["bbox"];
        af::ref<std::size_t> flags = temp["flags"];
        af::shared<std::size_t> selection;
        for (std::size_t i = 0; i < flags.size(); ++i) {
          flags[i] |= af::DontIntegrate;
          if (bbox[i][5] > frame[0]) {
            if (bbox[i][4] < frame[0]) {
              bbox[i][4] = frame[0];
            }
            selection.push_back(i);
          }
        }
        temp = select_rows_index(temp, selection.const_ref());
        extend(data, temp);
      }

      // Now select reflections from the block after
      if (index+1 < size()) {
        af::reflection_table temp = select_rows_index(
            lookup_.data(),
            lookup_.indices(index+1));
        af::ref<int6> bbox = temp["bbox"];
        af::ref<std::size_t> flags = temp["flags"];
        af::shared<std::size_t> selection;
        for (std::size_t i = 0; i < flags.size(); ++i) {
          flags[i] |= af::DontIntegrate;
          if (bbox[i][4] < frame[1]) {
            if (bbox[i][5] > frame[1]) {
              bbox[i][5] = frame[1];
            }
            selection.push_back(i);
          }
        }
        temp = select_rows_index(temp, selection.const_ref());
        extend(data, temp);
      }

      // Return the reflections
      return data;
    }

    /**
     * Accumulate the results.
     */
    void accumulate(std::size_t index, af::reflection_table result) {
      using namespace af::boost_python::flex_table_suite;
      DIALS_ASSERT(index < finished_.size());
      DIALS_ASSERT(finished_[index] == false);

      // Get the reflection indices
      af::const_ref<std::size_t> indices = lookup_.indices(index);

      // Get the total number of reflections to integrate
      std::size_t num_reflections = indices.size();

      // Check the number is less (because result includes adjacent reflections)
      DIALS_ASSERT(num_reflections < result.size());

      // Resize the input reflections to just those that were processed
      result.resize(indices.size());

      // Get the data array
      af::reflection_table data = lookup_.data();

      // Set the result in the lookup data
      set_selected_rows_index(data, indices, result);
      finished_[index] = true;
    }

  private:

    SimpleReflectionLookup lookup_;
    af::shared<bool> finished_;

  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H
