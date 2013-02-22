
#ifndef DIALS_SPOT_PREDICTION_BACKGROUND_INTENSITY_H
#define DIALS_SPOT_PREDICTION_BACKGROUND_INTENSITY_H

#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/math/mean_and_variance.h>
#include "../error.h"

namespace dials { namespace integration {

/** A class to calculate the background intensity of the reflection profile */
class BackgroundIntensity {

public:

    /** Default constructor */
    BackgroundIntensity()
        : delta_(0.1),
          max_iter_frac_(0.1) {}

    /**
     * Initialise the class with parameters
     * @param delta The deviation from normal
     * @param max_iter_frac The maximum numner of iterations as a fraction of
     *                      the elements in the input data array
     */
    BackgroundIntensity(double delta,
                        double max_iter_frac)
        : delta_(delta),
          max_iter_frac_(max_iter_frac)
    {
        DIALS_ASSERT(delta_ > 0.0 && delta_ < 1.0);
        DIALS_ASSERT(max_iter_frac_ > 0.0 && max_iter_frac_ < 1.0);
    }

    /*
     * Calculate the background intensity.
     *
     * Sort the pixels in order of ascending intensity. Then check if the
     * intensities are normally distributed. If not then remove the pixel
     * with the highest intensity from the list and check again. Keep going
     * untill the list of pixels is normally distributed, or the maximum
     * number of iterations is reached. Return the mean of the values as the
     * background intensity.
     *
     * This function modifies the input data
     *
     * @param pixels The list of pixels
     * @returns The background intensity value
     */
    double calculate(scitbx::af::flex_double_ref data)
    {
        // If no maximum is set, then set it to 0.1 N pixels if the max number
        // of iterations is greater than or equal the number of elements, then
        //
        int max_iter = (int)(max_iter_frac_ * data.size());

        // Sort the pixels into ascending intensity order
        std::sort(data.begin(), data.end());

        // Check if the data is normally distributed. If it is not, then remove
        // a value of high intensity and keep looping until it is. If the number
        // of iterations exceeds the maximum then exit the loop.
        int num_iter = 0;
        for (; num_iter < max_iter; ++num_iter) {
            if (this->is_data_normally_distributed(
                    scitbx::af::flex_double_const_ref(
                        data.begin(),
                        data.size()-num_iter))) {
                break;
            }
        }

        // Return the mean of the remaining pixels as the background intensity
        return scitbx::af::mean(scitbx::af::flex_double_const_ref(
                                    data.begin(),
                                    data.size()-num_iter));
    }

private:

    /**
     * Check if the data is normally distributed.
     *
     * Calculate the percentage number of data values within +- 1, 2 and 3
     * standard deviations of the mean value. Values that are normally
     * distributed will have:
     *      ~68.2% of points within +- 1 standard deviation
     *      ~95.4% of points within +- 2 standard deviation
     *      ~99.7% of points within +- 3 standard deviation
     *
     * If the percentages match (within +- delta) for the given data then the
     * function returns true, otherwise returns false.
     *
     * @param data The array of pixel values
     * @returns True/False
     */
    bool is_data_normally_distributed(const scitbx::af::flex_double_const_ref &data)
    {
        // Calculate the mean and standard deviation of the data
        int n_data = data.size();
        scitbx::math::mean_and_variance <double> mean_and_variance(
            scitbx::af::const_ref <double> (data.begin(), data.size()));
        double mean = mean_and_variance.mean();
        double sdev = mean_and_variance.unweighted_sample_standard_deviation();

        // Calculate the mean +- [1|2|3] standard deviations
        double sd1p = mean + 1 * sdev, sd1m = mean - 1 * sdev;
        double sd2p = mean + 2 * sdev, sd2m = mean - 2 * sdev;
        double sd3p = mean + 3 * sdev, sd3m = mean - 3 * sdev;

        // Calculate percentage of data within [1|2|3] standard deviations
        int n_perc1 = 0, n_perc2 = 0, n_perc3 = 0;
        for (int i = 0; i < n_data; ++i) {
            if (sd1m <= data[i] && data[i] < sd1p) n_perc1++;
            if (sd2m <= data[i] && data[i] < sd2p) n_perc2++;
            if (sd3m <= data[i] && data[i] < sd3p) n_perc3++;
        }
        double perc1 = (double)n_perc1 / (double)n_data;
        double perc2 = (double)n_perc2 / (double)n_data;
        double perc3 = (double)n_perc3 / (double)n_data;

        // Return whether the data is normally distributed
        return ((0.682 - delta_ <= perc1 && perc1 <= 0.682 + delta_) &&
                (0.954 - delta_ <= perc2 && perc2 <= 0.954 + delta_) &&
                (0.997 - delta_ <= perc3 && perc3 <= 0.997 + delta_));
    }

private:

    double delta_;
    double max_iter_frac_;
};

}} // namespace = dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_BACKGROUND_INTENSITY_H
