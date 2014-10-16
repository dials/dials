/*
 * reflection_table_statistics.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_REFLECTION_TABLE_STATISTICS_H
#define DIALS_ARRAY_FAMILY_REFLECTION_TABLE_STATISTICS_H

#include <dials/array_family/reflection_table.h>

namespace dials { namespace af {

  /**
   * A class to produce a summary of integation stats
   */
  class Summary {
  public:

    template <typename T>
    class MeanValue {
    public:
      MeanValue()
        : value_(0),
          count_(0) {}
      void operator+=(T v) {
        value_ += v;
        count_ ++;
      }
      operator T() const {
        if (count_ > 0) {
          return value_ / count_;
        }
        return 0;
      }
    private:
      T value_;
      std::size_t count_;
    };

    /**
     * A class to produce a summary of per image stats
     */
    class ImageSummary {
    public:

      class Element {
      public:
        int frame;
        std::size_t full;
        std::size_t part;
        MeanValue<double> sum_ios;
        MeanValue<double> prf_ios;
        Element()
          : frame(0),
            full(0),
            part(0) {}
        double sum_ios_value() const {
          return sum_ios;
        }
        double prf_ios_value() const {
          return prf_ios;
        }
      };

      class Data {
      public:
        Data() {}
        Data(int frame0, std::size_t nframes)
          : elements_(nframes) {
          for (std::size_t i = 0; i < elements_.size(); ++i) {
            elements_[i].frame = frame0 + i;
          }
        }
        Element& operator[](std::size_t index) {
          DIALS_ASSERT(index < elements_.size());
          return elements_[index];
        }
        Element at(std::size_t index) const {
          DIALS_ASSERT(index < elements_.size());
          return elements_[index];
        }
        std::size_t size() const {
          return elements_.size();
        }
      private:
        af::shared<Element> elements_;
      };

      ImageSummary(reflection_table self,
                   const af::const_ref<int2> &image_range) {

        // Check the input
        DIALS_ASSERT(image_range.size() > 0);
        DIALS_ASSERT(self.size() > 0);
        DIALS_ASSERT(self.is_consistent());
        DIALS_ASSERT(self.contains("bbox"));
        DIALS_ASSERT(self.contains("partiality"));
        DIALS_ASSERT(self.contains("intensity.sum.value"));
        DIALS_ASSERT(self.contains("intensity.sum.variance"));
        DIALS_ASSERT(self.contains("flags"));
        DIALS_ASSERT(self.contains("id"));

        // Get some data
        af::const_ref<int6> bbox = self["bbox"];
        af::const_ref<double> partiality = self["partiality"];
        af::const_ref<double> isum = self["intensity.sum.value"];
        af::const_ref<double> vsum = self["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = self["flags"];
        af::const_ref<std::size_t> id = self["id"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (self.contains("intensity.prf.value") &&
            self.contains("intensity.prf.variance")) {
          iprf = self["intensity.prf.value"];
          vprf = self["intensity.prf.variance"];
          prf = true;
        }

        // Allocate the arrays to the number of images
        for (std::size_t i = 0; i < image_range.size(); ++i) {
          DIALS_ASSERT(image_range[i][1] > image_range[i][0]);
          std::size_t nframes = image_range[i][1] - image_range[i][0];
          data_.push_back(Data(image_range[i][0], nframes));
        }

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the statistics
        for (std::size_t i = 0; i < self.size(); ++i) {

          // Get stuff for this experiment
          DIALS_ASSERT(id[i] < data_.size());
          int2 imrange = image_range[id[i]];
          Data& data = data_[id[i]];

          // Check the bbox is ok
          DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
          DIALS_ASSERT(bbox[i][4] >= imrange[0]);
          DIALS_ASSERT(bbox[i][5] <= imrange[1]);

          // Add to full or partial
          if (!(flags[i] & af::DontIntegrate)) {
            DIALS_ASSERT(partiality[i] >= 0 && partiality[i] <= 1);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - imrange[0];
              if (partiality[i] >= full_value) {
                data[k].full++;
              } else {
                data[k].part++;
              }
            }
          }

          // Add contribution for summation
          if (!(flags[i] & af::DontIntegrate) &&
               flags[i] & af::IntegratedSum &&
               vsum[i] > 0) {
            double ios = isum[i] / std::sqrt(vsum[i]);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - imrange[0];
              data[k].sum_ios += ios;
            }
          }

          // Add contribution for profile fitting
          if (!(flags[i] & af::DontIntegrate) &&
               flags[i] & af::IntegratedPrf &&
               prf && vprf[i] > 0) {
            double ios = iprf[i] / std::sqrt(vprf[i]);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - imrange[0];
              data[k].prf_ios += ios;
            }
          }
        }
      }

      Data data(std::size_t index) {
        return data_[index];
      }

      std::size_t size() const {
        return data_.size();
      }

    private:

      af::shared<Data> data_;
    };

    /**
     * Class to produce statistics at resolution bins
     */
    class ResolutionSummary {
    public:

      class Element {
      public:
        double2 d;
        MeanValue<double> sum_ios;
        MeanValue<double> sum_ios_full;
        MeanValue<double> sum_ios_part;
        MeanValue<double> prf_ios;
        MeanValue<double> prf_ios_full;
        MeanValue<double> prf_ios_part;
        Element() : d(0) {}
        double2 d_value() { return d; }
        double sum_ios_value() { return sum_ios; }
        double sum_ios_full_value() { return sum_ios_full; }
        double sum_ios_part_value() { return sum_ios_part; }
        double prf_ios_value() { return prf_ios; }
        double prf_ios_full_value() { return prf_ios_full; }
        double prf_ios_part_value() { return prf_ios_part; }
      };

      class Data {
      public:
        Data() {}
        Data(const af::const_ref<double> &d)
          : elements_(d.size()-1) {
          for (std::size_t i = 0; i < elements_.size(); ++i) {
            elements_[i].d = double2(d[i], d[i+1]);
          }
        }
        Element& operator[](std::size_t index) {
          DIALS_ASSERT(index < elements_.size());
          return elements_[index];
        }
        Element at(std::size_t index) const {
          DIALS_ASSERT(index < elements_.size());
          return elements_[index];
        }
        std::size_t size() const {
          return elements_.size();
        }
      private:
        af::shared<Element> elements_;
      };

      ResolutionSummary(reflection_table self,
                        std::size_t nbins,
                        std::size_t nexpr) {

        // Check the input
        DIALS_ASSERT(nbins > 0);
        DIALS_ASSERT(nexpr > 0);
        DIALS_ASSERT(self.size() > 0);
        DIALS_ASSERT(self.is_consistent());
        DIALS_ASSERT(self.contains("d"));
        DIALS_ASSERT(self.contains("partiality"));
        DIALS_ASSERT(self.contains("intensity.sum.value"));
        DIALS_ASSERT(self.contains("intensity.sum.variance"));
        DIALS_ASSERT(self.contains("flags"));
        DIALS_ASSERT(self.contains("id"));

        // Get some data
        af::const_ref<double> d = self["d"];
        af::const_ref<double> partiality = self["partiality"];
        af::const_ref<double> isum = self["intensity.sum.value"];
        af::const_ref<double> vsum = self["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = self["flags"];
        af::const_ref<std::size_t> id = self["id"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (self.contains("intensity.prf.value") &&
            self.contains("intensity.prf.variance")) {
          iprf = self["intensity.prf.value"];
          vprf = self["intensity.prf.variance"];
          prf = true;
        }

        // Compute the resolution bins
        double d_min = af::min(d);
        double d_max = af::max(d);
        DIALS_ASSERT(d_max > d_min);
        bin_range_ = (d_max - d_min) / nbins;
        for (std::size_t i = 0; i <= nbins; ++i) {
          bins_.push_back(d_min + bin_range_ * i);
        }

        // Allocate other arrays
        for (std::size_t i = 0; i < nexpr; ++i) {
          data_.push_back(Data(bins_.const_ref()));
        }

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the stats
        for (std::size_t i = 0; i < self.size(); ++i) {
          DIALS_ASSERT(id[i] < nexpr);
          Data& data = data_[id[i]];
          DIALS_ASSERT(d[i] >= d_min && d[i] <= d_max);
          int j = std::floor((d[i] - d_min) / bin_range_);
          DIALS_ASSERT(j >= 0 && j <= nbins);
          if (j == nbins) j = nbins-1;
          if (!(flags[i] & af::DontIntegrate)) {
            if (flags[i] & af::IntegratedSum &&
                vsum[i] > 0) {
              double ios = isum[i] / std::sqrt(vsum[i]);
              if (partiality[i] >= full_value) {
                data[j].sum_ios_full += ios;
              } else {
                data[j].sum_ios_part += ios;
              }
              data[j].sum_ios += ios;
            }
            if (flags[i] & af::IntegratedPrf && prf && vprf[i] > 0) {
              double ios = iprf[i] / std::sqrt(vprf[i]);
              if (partiality[i] >= full_value) {
                data[j].prf_ios_full += ios;
              } else {
                data[j].prf_ios_part += ios;
              }
              data[j].prf_ios += ios;
            }
          }
        }
      }

      Data data(std::size_t index) const {
        return data_[index];
      }

      std::size_t size() const {
        return data_.size();
      }

    private:

      double bin_range_;
      af::shared<double> bins_;
      af::shared<Data> data_;
    };

    /**
     * A class to produce summary stats for whole dataset
     */
    class WholeSummary {
    public:


      class Data {
      public:

        MeanValue<double> sum_ios;
        MeanValue<double> sum_ios_full;
        MeanValue<double> sum_ios_part;
        MeanValue<double> prf_ios;
        MeanValue<double> prf_ios_full;
        MeanValue<double> prf_ios_part;
        double sum_ios_value() { return sum_ios; }
        double sum_ios_full_value() { return sum_ios_full; }
        double sum_ios_part_value() { return sum_ios_part; }
        double prf_ios_value() { return prf_ios; }
        double prf_ios_full_value() { return prf_ios_full; }
        double prf_ios_part_value() { return prf_ios_part; }
      };

      WholeSummary(reflection_table self, std::size_t nexpr)
        : data_(nexpr) {

        // Check the input
        DIALS_ASSERT(nexpr > 0);
        DIALS_ASSERT(self.size() > 0);
        DIALS_ASSERT(self.is_consistent());
        DIALS_ASSERT(self.contains("partiality"));
        DIALS_ASSERT(self.contains("intensity.sum.value"));
        DIALS_ASSERT(self.contains("intensity.sum.variance"));
        DIALS_ASSERT(self.contains("flags"));
        DIALS_ASSERT(self.contains("id"));

        // Get some data
        af::const_ref<double> partiality = self["partiality"];
        af::const_ref<double> isum = self["intensity.sum.value"];
        af::const_ref<double> vsum = self["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = self["flags"];
        af::const_ref<std::size_t> id = self["id"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (self.contains("intensity.prf.value") &&
            self.contains("intensity.prf.variance")) {
          iprf = self["intensity.prf.value"];
          vprf = self["intensity.prf.variance"];
          prf = true;
        }

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the stats
        for (std::size_t i = 0; i < self.size(); ++i) {
          std::size_t j = id[i];
          DIALS_ASSERT(j < nexpr);
          if (!(flags[i] & af::DontIntegrate)) {
            if (flags[i] & af::IntegratedSum && vsum[i] > 0) {
              double ios = isum[i] / std::sqrt(vsum[i]);
              if (partiality[i] >= full_value) {
                data_[j].sum_ios_full += ios;
              } else {
                data_[j].sum_ios_part += ios;
              }
              data_[j].sum_ios += ios;
            }
            if (flags[i] & af::IntegratedPrf && prf && vprf[i] > 0) {
              double ios = iprf[i] / std::sqrt(vprf[i]);
              if (partiality[i] >= full_value) {
                data_[j].prf_ios_full += ios;
              } else {
                data_[j].prf_ios_part += ios;
              }
              data_[j].prf_ios += ios;
            }
          }
        }
      }

      Data data(std::size_t index) const {
        DIALS_ASSERT(index < data_.size());
        return data_[index];
      }

      std::size_t size() const {
        return data_.size();
      }

    private:
      af::shared<Data> data_;
    };

    Summary(reflection_table data,
            const af::const_ref<int2> &image_range,
            std::size_t nresbins)
       :img_summary_(data, image_range),
        res_summary_(data, nresbins, image_range.size()),
        who_summary_(data, image_range.size()) {
    }

    ImageSummary image_summary() {
      return img_summary_;
    }

    ResolutionSummary resolution_summary() {
      return res_summary_;
    }

    WholeSummary whole_summary() {
      return who_summary_;
    }

  private:

    ImageSummary img_summary_;
    ResolutionSummary res_summary_;
    WholeSummary who_summary_;
  };

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_REFLECTION_TABLE_STATISTICS_H
