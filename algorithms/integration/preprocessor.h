/*
 * preprocessor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PREPROCESSOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PREPROCESSOR_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctbx/miller/index_generator.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>


namespace dials { namespace algorithms {

  /**
   * A class to do powder ring filtering.
   */
  class PowderRingFilter {
  public:

    /**
     * @param uc The unit cell
     * @param sg The space group
     * @param d_min The maximum resolution to filter at
     * @param d_max The minimum resolution to filter at
     * @param width The width of the filter
     */
    PowderRingFilter(
        cctbx::uctbx::unit_cell uc,
        cctbx::sgtbx::space_group sg,
        double d_min,
        double width)
      : d_min_(d_min),
        width_(width),
        unit_cell_(uc),
        space_group_(sg) {
      DIALS_ASSERT(d_min > 0);
      DIALS_ASSERT(width > 0);

      // Correct unit cell
      uc = sg.average_unit_cell(uc);

      // Generate a load of indices
      cctbx::miller::index_generator generator(uc, sg.type(), false, d_min);
      af::shared< cctbx::miller::index<> > indices = generator.to_array();

      // Calculate the d spacings
      for (std::size_t i = 0; i < indices.size(); ++i) {
        d_spacings_.push_back(uc.d(indices[i]));
      }

      // Sort the d spacings by resolution
      std::sort(d_spacings_.begin(), d_spacings_.end());
    }

    /**
     * @returns True if within powder ring.
     */
    bool operator()(double d) const {
      double half_width = width_ / 2.0;
      for (std::size_t i = 0; i < d_spacings_.size(); ++i) {
        if (std::abs(d - d_spacings_[i]) < half_width) {
          return true;
        }
      }
      return false;
    }

    /** @returns The d spacings used */
    af::shared<double> d_spacings() const {
      return d_spacings_;
    }

    /** @returns The width used */
    double width() const {
      return width_;
    }

    /** @returns The maximum resolution to filter at */
    double d_min() const {
      return d_min_;
    }

    /** @returns The unit cell */
    cctbx::uctbx::unit_cell unit_cell() const {
      return unit_cell_;
    }

    /** @returns The space group */
    cctbx::sgtbx::space_group space_group() const {
      return space_group_;
    }

  private:

    double d_min_;
    double d_max_;
    double width_;
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group space_group_;
    af::shared<double> d_spacings_;
  };


  /**
   * A class to encapsulate multiple powder ring filters
   */
  class MultiPowderRingFilter {
  public:

    MultiPowderRingFilter() {}

    /**
     * Add another powder ring filter.
     */
    void add(const PowderRingFilter &filter) {
      filters_.push_back(filter);
    }

    /**
     * @returns The powder ring filter at index
     */
    const PowderRingFilter& operator[](std::size_t index) const {
      DIALS_ASSERT(index < filters_.size());
      return filters_[index];
    }

    /**
     * @returns The powder ring filter at index
     */
    PowderRingFilter& operator[](std::size_t index) {
      DIALS_ASSERT(index < filters_.size());
      return filters_[index];
    }

    /**
     * @returns True if within powder ring.
     */
    bool operator()(double d) const {
      for (std::size_t i = 0; i < filters_.size(); ++i) {
        if (filters_[i](d)) {
          return true;
        }
      }
      return false;
    }

    /** @returns The number of filters */
    std::size_t size() const {
      return filters_.size();
    }

  private:
    af::shared<PowderRingFilter> filters_;
  };


  /**
   * Class to preprocess reflections for integration.
   */
  class Preprocessor {
  public:

    /**
     * Setup and perform the preprocessing.
     * @param reflections The reflections to process.
     * @param filter The powder ring filter
     * @param min_zeta The minimum allowed zeta
     */
    Preprocessor(
            af::reflection_table reflections,
            const MultiPowderRingFilter &filter,
            double min_zeta)
        : num_total_(0),
          num_strong_(0),
          num_filtered_(0),
          num_integrate_(0),
          num_overlap_bg_(0),
          num_overlap_fg_(0),
          num_inpowder_(0),
          nums_overlap_bg_(0),
          nums_overlap_fg_(0),
          nums_inpowder_(0),
          nums_reference_(0),
          small_x_(0),
          small_y_(0),
          small_z_(0),
          large_x_(0),
          large_y_(0),
          large_z_(0),
          filter_(filter),
          min_zeta_(min_zeta) {

      // Check the input
      DIALS_ASSERT(min_zeta > 0);

      // Check reflection table contains expected properties
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.contains("flags"));
      DIALS_ASSERT(reflections.contains("id"));
      DIALS_ASSERT(reflections.contains("panel"));
      DIALS_ASSERT(reflections.contains("miller_index"));
      DIALS_ASSERT(reflections.contains("entering"));
      DIALS_ASSERT(reflections.contains("s1"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(reflections.contains("xyzcal.mm"));
      DIALS_ASSERT(reflections.contains("bbox"));
      DIALS_ASSERT(reflections.contains("zeta"));
      DIALS_ASSERT(reflections.contains("d"));

      // Get some arrays needed for preprocessing
      af::ref<std::size_t> flags = reflections["flags"];
      af::const_ref<int6> bbox = reflections["bbox"];
      af::const_ref<double> zeta = reflections["zeta"];
      af::const_ref<double> d = reflections["d"];

      // Compute some summary stuff
      num_total_ = reflections.size();
      for (std::size_t i = 0; i < bbox.size(); ++i) {

        // Reset any filtering and decision flags
        flags[i] &= ~af::DontIntegrate;
        flags[i] &= ~af::InPowderRing;
        flags[i] &= ~af::ReferenceSpot;

        // Filter the reflections by zeta
        if (std::abs(zeta[i]) < min_zeta) {
          flags[i] |= af::DontIntegrate;
        }

        // Check if reflections are predicted on powder rings
        if (filter(d[i])) {
          flags[i] |= af::InPowderRing;
        }

        // If a spot is strong and neither of the bad flags are set
        // then set the strong spot as a spot to use for reference profiles
        if (flags[i] & af::Strong) {
          std::size_t code = af::DontIntegrate | af::InPowderRing;
          if (!(flags[i] | code)) {
            flags[i] |= af::ReferenceSpot;
          }
        }

        // Update some counters
        num_strong_ += flags[i] & af::Strong ? 1 : 0;
        num_filtered_ += flags[i] & af::DontIntegrate ? 1 : 0;
        num_inpowder_ += flags[i] & af::InPowderRing ? 1 : 0;
        if (flags[i] & af::Strong) {
          nums_inpowder_ += flags[i] & af::InPowderRing ? 1 : 0;
          nums_reference_ += flags[i] & af::ReferenceSpot ? 1 : 0;
        } else {
          DIALS_ASSERT(!(flags[i] & af::ReferenceSpot));
        }

        // Compute stuff about the bounding box
        int6 b = bbox[i];
        DIALS_ASSERT(b[1] > b[0]);
        DIALS_ASSERT(b[3] > b[2]);
        DIALS_ASSERT(b[5] > b[4]);
        if (i == 0) {
          small_x_ = large_x_ = (std::size_t)(b[1] - b[0]);
          small_y_ = large_y_ = (std::size_t)(b[3] - b[2]);
          small_z_ = large_z_ = (std::size_t)(b[5] - b[4]);
        } else {
          small_x_ = std::min(small_x_, (std::size_t)(b[1] - b[0]));
          small_y_ = std::min(small_y_, (std::size_t)(b[3] - b[2]));
          small_z_ = std::min(small_z_, (std::size_t)(b[5] - b[4]));
          large_x_ = std::max(large_x_, (std::size_t)(b[1] - b[0]));
          large_y_ = std::max(large_y_, (std::size_t)(b[3] - b[2]));
          large_z_ = std::max(large_z_, (std::size_t)(b[5] - b[4]));
        }
      }
      num_integrate_ = num_total_ - num_filtered_;
    }

    /**
     * @returns a text summary of the preprocessing
     */
    std::string summary() const {

      // Helper to print a percent
      struct percent {
        double p_;
        percent(double p) : p_(p) {}
        std::string s() const {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(2) << p_ << "%";
          return ss.str();
        }
      };

      // Helper to print a number and percent
      struct number_and_percent {
        std::size_t t_;
        number_and_percent(std::size_t t) : t_(t) {}
        std::string s(std::size_t n) const {
          DIALS_ASSERT(n >= 0 && n <= t_);
          double p = t_ > 0 ? (double)n / (double)t_ : 0.0;
          std::stringstream ss;
          ss << n << " (" << percent(100.0 * p).s() << ") ";
          return ss.str();
        }
      };

      // Format the summary output
      number_and_percent nt(num_total_);
      number_and_percent ns(num_strong_);
      std::stringstream ss;
      ss
        << "Preprocessed reflections:\n"
        << "\n"
        << filter_info_string()
        << " Filtered reflections with min zeta =    " << min_zeta_ << "\n"
        << "\n"
        << " Number of reflections:                  " << num_total_ << "\n"
        << " Number of strong reflections:           " << num_strong_ << "\n"
        << " Number filtered with zeta:              " << nt.s(num_filtered_) << "\n"
        << " Number of reflections to integrate:     " << nt.s(num_integrate_) << "\n"
        << " Number recorded on powder rings:        " << nt.s(num_inpowder_) << "\n"
        << " Number strong recorded on powder rings: " << ns.s(nums_inpowder_) << "\n"
        << " Number strong for reference creation:   " << ns.s(nums_reference_) << "\n"
        << " Smallest measurement box size (x):      " << small_x_ << "\n"
        << " Smallest measurement box size (y):      " << small_y_ << "\n"
        << " Smallest measurement box size (z):      " << small_z_ << "\n"
        << " Largest measurement box size (x):       " << large_x_ << "\n"
        << " Largest measurement box size (y):       " << large_y_ << "\n"
        << " Largest measurement box size (z):       " << large_z_ << "\n"
        ;
      return ss.str();
    }

  private:

    /** @returns A summary of the input filters */
    std::string filter_info_string() const {
      std::stringstream ss_out;
      for (std::size_t i = 0; i < filter_.size(); ++i) {
        const PowderRingFilter &f = filter_[i];
        std::stringstream ss;
        ss << " Identified reflections on the following powder rings";
        ss << "\n";
        ss << "  Unit cell:   ";
        af::double6 p = f.unit_cell().parameters();
        ss << p[0];
        for (std::size_t j = 1; j < p.size(); ++j) {
          ss << ", " << p[j];
        }
        ss << "\n";
        ss << "  Space group: " << f.space_group().type().number() << "\n";
        ss << "  D min:       " << f.d_min() << "\n";
        ss << "  D width:     " << f.width() << "\n";
        ss << "  D spacings:  ";
        af::shared<double> d = f.d_spacings();
        for (std::size_t j = 0; j < d.size(); ++j) {
          if (j % 6 == 0) {
            ss << "\n";
            ss << "   ";
          }
          ss << std::fixed << std::setprecision(5) << d[j] << " ";
        }
        ss << "\n";
        ss << "\n";
        ss_out << ss.str();
      }
      return ss_out.str();
    }

    std::size_t num_total_;
    std::size_t num_strong_;
    std::size_t num_filtered_;
    std::size_t num_integrate_;
    std::size_t num_overlap_bg_;
    std::size_t num_overlap_fg_;
    std::size_t num_inpowder_;
    std::size_t nums_overlap_bg_;
    std::size_t nums_overlap_fg_;
    std::size_t nums_inpowder_;
    std::size_t nums_reference_;
    std::size_t small_x_;
    std::size_t small_y_;
    std::size_t small_z_;
    std::size_t large_x_;
    std::size_t large_y_;
    std::size_t large_z_;
    MultiPowderRingFilter filter_;
    double min_zeta_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PREPROCESSOR_H
