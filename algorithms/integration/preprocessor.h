
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
        width_(width) {
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

  private:

    double d_min_;
    double d_max_;
    double width_;
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


  class Preprocessor {
  public:

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
          num_icering_(0),
          nums_overlap_bg_(0),
          nums_overlap_fg_(0),
          nums_icering_(0),
          nums_reference_(0),
          small_x_(0),
          small_y_(0),
          small_z_(0),
          large_x_(0),
          large_y_(0),
          large_z_(0) {

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
      af::const_ref<std::size_t> id = reflections["id"];
      af::const_ref< vec3<double> > s1 = reflections["s1"];
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
        if (std::abs(zeta[i]) < 0.05) {
          flags[i] |= af::DontIntegrate;
        }

        // Check if reflections are predicted on ice rings

        // Check the flags
        if (flags[i] & af::Strong) {
          num_strong_++;
        }
        if (flags[i] & af::DontIntegrate) {
          num_filtered_++;
        } else {
          num_integrate_++;
        }
        if (flags[i] & af::ReferenceSpot) {
          nums_reference_++;
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
    }

    std::string summary() const {

      // Helper to print a percent
      struct percent {
        double p_;
        percent(double p) : p_(p) {}
        std::string s() const {
          std::stringstream ss;
          ss << std::setprecision(2) << p_ << "%";
          return ss.str();
        }
      };

      // Helper to print a number and percent
      struct number_and_percent {
        std::size_t t_;
        number_and_percent(std::size_t t) : t_(t) {}
        std::string s(std::size_t n) const {
          DIALS_ASSERT(n >= 0 && n <= t_);
          double p = t_ > 0 ? n / t_ : 0.0;
          std::stringstream ss;
          ss << n << " (" << percent(p).s() << ") ";
          return ss.str();
        }
      };

      // Format the summary output
      number_and_percent nt(num_total_);
      number_and_percent ns(num_strong_);
      std::stringstream ss;
      ss <<
        "Preprocessed reflections:\n"
        "\n"
        " Number of reflections:                  " << num_total_ << "\n"
        " Number of strong reflections:           " << num_strong_ << "\n"
        " Number filtered with zeta:              " << nt.s(num_filtered_) << "\n"
        " Number of reflection to integrate:      " << nt.s(num_integrate_) << "\n"
        " Number overlapping (background):        " << nt.s(num_overlap_bg_) << "\n"
        " Number overlapping (foreground):        " << nt.s(num_overlap_fg_) << "\n"
        " Number recorded on ice rings:           " << nt.s(num_icering_) << "\n"
        " Number strong overlapping (background): " << ns.s(nums_overlap_bg_) << "\n"
        " Number strong overlapping (foreground): " << ns.s(nums_overlap_fg_) << "\n"
        " Number strong recorded on ice rings:    " << ns.s(nums_icering_) << "\n"
        " Number strong for reference creation:   " << ns.s(nums_reference_) << "\n"
        " Smallest measurement box size (x):      " << small_x_ << "\n"
        " Smallest measurement box size (y):      " << small_y_ << "\n"
        " Smallest measurement box size (z):      " << small_z_ << "\n"
        " Largest measurement box size (x):       " << large_x_ << "\n"
        " Largest measurement box size (y):       " << large_y_ << "\n"
        " Largest measurement box size (z):       " << large_z_ << "\n"
        ;
      return ss.str();
    }

  private:
    std::size_t num_total_;
    std::size_t num_strong_;
    std::size_t num_filtered_;
    std::size_t num_integrate_;
    std::size_t num_overlap_bg_;
    std::size_t num_overlap_fg_;
    std::size_t num_icering_;
    std::size_t nums_overlap_bg_;
    std::size_t nums_overlap_fg_;
    std::size_t nums_icering_;
    std::size_t nums_reference_;
    std::size_t small_x_;
    std::size_t small_y_;
    std::size_t small_z_;
    std::size_t large_x_;
    std::size_t large_y_;
    std::size_t large_z_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PREPROCESSOR_H
