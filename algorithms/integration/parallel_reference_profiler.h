/*
 * parallel_parallel_reference_profiler.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_PARALLEL_REFERENCE_PROFILER_H
#define DIALS_ALGORITHMS_INTEGRATION_PARALLEL_REFERENCE_PROFILER_H

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
#include <dials/algorithms/integration/parallel_integrator.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/centroid/centroid.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <map>

#include <dials/algorithms/integration/interfaces.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Panel;
  using dxtbx::model::Scan;

  using dxtbx::ImageSweep;
  using dxtbx::format::Image;
  using dxtbx::format::ImageTile;

  using dials::model::AdjacencyList;
  using dials::model::Shoebox;

  /**
   * A class to integrate a single reflection
   */
  class ReflectionReferenceProfiler {
  public:
    /**
     * Initialise the parallel_reference_profiler
     * @param compute_mask The mask calculation function
     * @param compute_background The background calculation function
     * @param compute_reference The intensity calculation function
     * @param buffer The image buffer array
     * @param zstart The first image index
     * @param underload The underload value
     * @param overload The overload value
     */
    ReflectionReferenceProfiler(const MaskCalculatorIface &compute_mask,
                                const BackgroundCalculatorIface &compute_background,
                                ReferenceCalculatorIface &compute_reference,
                                const Buffer &buffer,
                                int zstart,
                                double underload,
                                double overload,
                                bool debug)
        : compute_mask_(compute_mask),
          compute_background_(compute_background),
          compute_reference_(compute_reference),
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
    void operator()(std::size_t index,
                    af::ref<af::Reflection> reflection_list,
                    const AdjacencyList &adjacency_list) const {
      af::Reflection reflection;
      std::vector<af::Reflection> adjacent_reflections;

      // Get the reflection data
      get_reflection(
        index, reflection_list, adjacency_list, reflection, adjacent_reflections);

      // Extract the shoebox data
      extract_shoebox(buffer_, reflection, zstart_, underload_, overload_);

      // Compute the mask
      compute_mask_(reflection);

      // Set all the bounding boxes of adjacent reflections
      // And compute the mask for these reflections too.
      for (std::size_t i = 0; i < adjacent_reflections.size(); ++i) {
        adjacent_reflections[i]["bbox"] = reflection.get<int6>("bbox");
        adjacent_reflections[i]["shoebox"] = reflection.get<Shoebox<> >("shoebox");
        compute_mask_(adjacent_reflections[i], true);
      }

      // Compute the background
      try {
        compute_background_(reflection);
      } catch (dials::error) {
        finalize_shoebox(reflection, adjacent_reflections, underload_, overload_);
        return;
      }

      // Compute the centroid
      compute_centroid(reflection);

      // Compute the summed intensity
      compute_summed_intensity(reflection);

      // Compute the profile fitted intensity
      try {
        compute_reference_(reflection);
      } catch (dials::error) {
        // pass
      }

      // Erase the shoebox
      finalize_shoebox(reflection, adjacent_reflections, underload_, overload_);

      // Set the reflection data
      set_reflection(index, reflection_list, reflection);
    }

  protected:
    /**
     * Get the reflection data in a thread safe manner
     * @param index The reflection index
     * @param reflection_list The reflection list
     * @param adjacency_list The adjacency list
     * @param reflection The reflection data
     * @param adjacent_reflections The adjacent reflections
     */
    void get_reflection(std::size_t index,
                        const af::const_ref<af::Reflection> &reflection_list,
                        const AdjacencyList &adjacency_list,
                        af::Reflection &reflection,
                        std::vector<af::Reflection> &adjacent_reflections) const {
      DIALS_ASSERT(index < reflection_list.size());

      // Get the lock
      boost::lock_guard<boost::mutex> guard(mutex_);

      // Get the reflection
      reflection = reflection_list[index];

      // Get the adjacent reflections
      adjacent_reflections.reserve(adjacency_list.vertex_num_edges(index));
      AdjacencyList::edge_iterator_range edges = adjacency_list.edges(index);
      for (AdjacencyList::edge_iterator it = edges.first; it != edges.second; ++it) {
        DIALS_ASSERT(it->first == index);
        DIALS_ASSERT(it->second < reflection_list.size());
        adjacent_reflections.push_back(reflection_list[it->second]);
      }
    }

    /**
     * Set the reflection data in a thread safe way
     * @param index The reflection index
     * @param reflection_list The reflection list
     * @param reflection The reflection data
     */
    void set_reflection(std::size_t index,
                        af::ref<af::Reflection> reflection_list,
                        const af::Reflection &reflection) const {
      DIALS_ASSERT(index < reflection_list.size());
      boost::lock_guard<boost::mutex> guard(mutex_);
      reflection_list[index] = reflection;
    }

    /**
     * Before exiting do some stuff on the shoebox
     * @param reflection The reflection
     * @param adjacent_reflections The adjancent reflections
     */
    void finalize_shoebox(af::Reflection &reflection,
                          std::vector<af::Reflection> &adjacent_reflections,
                          double underload,
                          double overload) const {
      // Inspect the pixels
      inspect_pixels(reflection, underload, overload);

      // Delete the shoebox
      delete_shoebox(reflection, adjacent_reflections);
    }

    /**
     * Inspect the pixel and mask values
     * @param reflection The reflection
     */
    void inspect_pixels(af::Reflection &reflection,
                        double underload,
                        double overload) const {
      typedef Shoebox<>::float_type float_type;

      // Get the shoebox
      Shoebox<> &sbox = reflection.get<Shoebox<> >("shoebox");
      std::size_t flags = reflection.get<std::size_t>("flags");

      // Get the pixel data
      af::const_ref<float_type, af::c_grid<3> > data = sbox.data.const_ref();
      af::const_ref<int, af::c_grid<3> > mask = sbox.mask.const_ref();
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

      // Check pixel values
      std::size_t n_valid = 0;
      std::size_t n_background = 0;
      std::size_t n_background_used = 0;
      std::size_t n_foreground = 0;
      std::size_t mask_code1 = Valid;
      std::size_t mask_code2 = Valid | Background;
      std::size_t mask_code3 = Valid | Background | BackgroundUsed;
      std::size_t mask_code4 = Valid | Foreground;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        double d = data[i];
        int m = mask[i];

        if (d >= overload) {
          flags |= af::Overloaded;
        }

        if ((m & Background) && !(m & Valid)) {
          flags |= af::BackgroundIncludesBadPixels;
        }

        if ((m & Foreground) && !(m & Valid)) {
          flags |= af::ForegroundIncludesBadPixels;
        }

        if ((m & Background) && (m & Overlapped)) {
          flags |= af::OverlappedBg;
        }

        if ((m & Foreground) && (m & Overlapped)) {
          flags |= af::OverlappedFg;
        }

        if ((m & mask_code1) == mask_code1) {
          n_valid++;
        }

        if ((m & mask_code2) == mask_code2) {
          n_background++;
        }

        if ((m & mask_code3) == mask_code3) {
          n_background_used++;
        }

        if ((m & mask_code4) == mask_code4) {
          n_foreground++;
        }
      }

      // Set some information in the reflection
      reflection["num_pixels.valid"] = (int)n_valid;
      reflection["num_pixels.background"] = (int)n_background;
      reflection["num_pixels.background_used"] = (int)n_background_used;
      reflection["num_pixels.foreground"] = (int)n_foreground;
      reflection["flags"] = flags;
    }

    /**
     * Delete the shoebox
     * @param reflection The reflection
     * @param adjacent_reflections The adjancent reflections
     */
    void delete_shoebox(af::Reflection &reflection,
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
    void extract_shoebox(const Buffer &buffer,
                         af::Reflection &reflection,
                         int zstart,
                         double underload,
                         double overload) const {
      typedef af::const_ref<Buffer::float_type, af::c_grid<2> > data_buffer_type;
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref<float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref<int, af::c_grid<3> > mask = shoebox.mask.ref();
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
        int kk = z0 + k - zstart;
        if (kk < 0 || kk >= buffer.num_images()) {
          continue;
        }
        data_buffer_type data_buffer = buffer.data(panel, kk);
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int jj = y0 + j;
            int ii = x0 + i;
            if (jj >= 0 && ii >= 0 && jj < data_buffer.accessor()[0]
                && ii < data_buffer.accessor()[1]) {
              double d = data_buffer(jj, ii);
              int m = (d > underload && d < overload) ? Valid : 0;
              data(k, j, i) = d;
              mask(k, j, i) = m;
            } else {
              data(k, j, i) = 0;
              mask(k, j, i) = 0;
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
      Shoebox<> shoebox = reflection.get<Shoebox<> >("shoebox");
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
      Shoebox<> shoebox = reflection.get<Shoebox<> >("shoebox");
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
    ReferenceCalculatorIface &compute_reference_;
    const Buffer &buffer_;
    int zstart_;
    double underload_;
    double overload_;
    bool debug_;
    mutable boost::mutex mutex_;
  };

  /**
   * A class to perform parallel integration
   */
  class ParallelReferenceProfiler {
  public:
    /**
     * Do the integration
     * @param reflections The reflection table
     * @param imageset The imageset
     * @param compute_mask The mask calulcation function
     * @param compute_background The background calculation function
     * @param compute_reference The intensity calculation function
     * @param logger The logger class
     * @param nthreads The number of parallel threads
     * @param buffer_size The buffer_size
     * @param use_dynamic_mask Use the dynamic mask if present
     * @param debug Add debug output
     */
    ParallelReferenceProfiler(af::reflection_table reflections,
                              ImageSweep imageset,
                              const MaskCalculatorIface &compute_mask,
                              const BackgroundCalculatorIface &compute_background,
                              ReferenceCalculatorIface &compute_reference,
                              const Logger &logger,
                              std::size_t nthreads,
                              std::size_t buffer_size,
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
      if (buffer_size == 0 || buffer_size > zsize) {
        buffer_size = zsize;
      }

      // Get the starting frame and the underload/overload values
      int zstart = scan.get_array_range()[0];
      double underload = detector[0].get_trusted_range()[0];
      double overload = detector[0].get_trusted_range()[1];
      DIALS_ASSERT(underload < overload);
      for (std::size_t i = 1; i < detector.size(); ++i) {
        DIALS_ASSERT(underload == detector[i].get_trusted_range()[0]);
        DIALS_ASSERT(overload == detector[i].get_trusted_range()[1]);
      }

      // Get the reflection flags and bbox
      af::const_ref<std::size_t> panel =
        reflections.get<std::size_t>("panel").const_ref();
      af::const_ref<int6> bbox = reflections.get<int6>("bbox").const_ref();
      af::ref<std::size_t> flags = reflections.get<std::size_t>("flags").ref();

      // Reset the flags
      reset_flags(flags);

      // Find the overlapping reflections
      AdjacencyList overlaps = find_overlapping_multi_panel(bbox, panel);

      // Allocate the array for the image data
      Buffer buffer(
        detector, zsize, buffer_size, underload, imageset.get_static_mask());

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

      // Create the reflection parallel_reference_profiler. This class is called for
      // each reflection to integrate the data
      ReflectionReferenceProfiler parallel_reference_profiler(compute_mask,
                                                              compute_background,
                                                              compute_reference,
                                                              buffer,
                                                              zstart,
                                                              underload,
                                                              overload,
                                                              debug);

      // Do the integration
      process(lookup,
              parallel_reference_profiler,
              buffer,
              reflection_array.ref(),
              overlaps,
              imageset,
              bbox,
              flags,
              nthreads,
              use_dynamic_mask,
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
     */
    static std::size_t compute_required_memory(ImageSweep imageset,
                                               std::size_t block_size) {
      DIALS_ASSERT(imageset.get_detector() != NULL);
      DIALS_ASSERT(imageset.get_scan() != NULL);
      Detector detector = *imageset.get_detector();
      Scan scan = *imageset.get_scan();
      block_size = std::min(block_size, (std::size_t)scan.get_num_images());
      std::size_t nelements = 0;
      for (std::size_t i = 0; i < detector.size(); ++i) {
        std::size_t xsize = detector[i].get_image_size()[0];
        std::size_t ysize = detector[i].get_image_size()[1];
        nelements += xsize * ysize;
      }
      nelements *= block_size;
      std::size_t nbytes = nelements * sizeof(double);
      return nbytes;
    }

    /**
     * Static method to get the memory in bytes needed
     * @param imageset the imageset class
     * @param max_memory_usage The maximum memory usage
     */
    static std::size_t compute_max_block_size(ImageSweep imageset,
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
        flags[i] &= ~af::UsedInModelling;
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
    void process(const Lookup &lookup,
                 const ReflectionReferenceProfiler &parallel_reference_profiler,
                 Buffer &buffer,
                 af::ref<af::Reflection> reflections,
                 const AdjacencyList &overlaps,
                 ImageSweep imageset,
                 af::const_ref<int6> bbox,
                 af::const_ref<std::size_t> flags,
                 std::size_t nthreads,
                 bool use_dynamic_mask,
                 const Logger &logger) const {
      using dials::util::ThreadPool;

      // Create the thread pool
      ThreadPool pool(nthreads);

      // Get the size of the array
      int zstart = imageset.get_scan()->get_array_range()[0];
      std::size_t zsize = imageset.size();

      // Create the buffer manager
      BufferManager bm(buffer, bbox, flags, zstart);

      // Loop through all the images
      for (std::size_t i = 0; i < zsize; ++i) {
        // Copy the image to the buffer. If the image number is greater than the
        // buffer size (i.e. we are now deleting old images) then wait for the
        // threads to finish so that we don't end up reading the wrong data
        if (imageset.is_marked_for_rejection(i)) {
          bm.copy_when_ready(imageset.get_corrected_data(i), false, i);
        } else if (use_dynamic_mask) {
          bm.copy_when_ready(
            imageset.get_corrected_data(i), imageset.get_dynamic_mask(i), i);
        } else {
          bm.copy_when_ready(imageset.get_corrected_data(i), i);
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
          bm.post(
            pool,
            boost::bind(&ReflectionReferenceProfiler::operator(),
                        boost::ref(parallel_reference_profiler),
                        k,
                        af::ref<af::Reflection>(&reflections[0], reflections.size()),
                        boost::ref(overlaps)),
            bbox[k][4]);
        }

        // Print some output
        std::ostringstream ss;
        ss << "Modelling " << std::setw(5) << count << " reflections on image "
           << std::setw(6) << zstart + i;
        logger.info(ss.str().c_str());
      }

      // Wait for all the integration jobs to complete
      bm.wait(pool);
    }

    af::reflection_table reflections_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_PARALLEL_REFERENCE_PROFILER_H
