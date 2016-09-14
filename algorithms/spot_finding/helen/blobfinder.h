/*
 * blobfinder.h
 *
 * An exploratory coding adventure into finding spot-shaped blobs on
 * XFEL detectors.
 *
 * Copyright © 2016 Helen Ginn
 *
 * Authors:
 *   2016           Helen Ginn <helen@strubi.ox.ac.uk>
 *
 * Notes:
 *
 * Nobody's told me that I can't put opening braces on the next line
 * so you're just going to have to deal with it. Haha!
 */

/* ------------------------------------------------------------------------
 * Definitions
 * ------------------------------------------------------------------------*/

#ifndef BLOBFINDER_H
#define BLOBFINDER_H

/* Dummy intensity to scale the 'ideal spot' shape */
#define BLOB_MAX_INTENSITY 300

/* To check the code is working, print out various bits and pieces */
#define BLOB_DEBUG 1

/* When sanity checking a potential spot, what is the maximum proportion
 * of masked pixels we are willing to accept?
 **/
#define BLOB_MAX_MASKED_PIXELS 0.5

#ifdef BLOB_DEBUG
#include <iostream>
#endif

#include <cmath>

/* ------------------------------------------------------------------------
 * functions called under the core functions, still specialised (Level 3)
 * ------------------------------------------------------------------------*/


class BlobThresholdAlgorithm
{
private:
        int _pixels_per_row;
        int _row_count;
        int _image_pix_num;
        int _exp_spot_dimension;
        int _model_half_dim;
        double _min_blob_score;

        int *_ideal_spot;
        int *_shifts;
        int _spot_pix_num;

        /* for now, hoping to remove this parameter */
        double _global_threshold;

        /**
         * Calculates what an 'ideal' spot looks like and stores this in
         * the member _ideal_spot[].
         * Also calculates the 'shifts' from 0 assuming the centre of an
         * odd-pixel diameter spot is in the centre, and for an even-diameter
         * spot, displaced to the bottom right slightly.
         **/
        void calculate_ideal_spot()
        {
                /* add two "background" pixels for comparison */
                int model_diameter = _exp_spot_dimension + 2;

                _spot_pix_num = pow(model_diameter, 2);

                _ideal_spot = new int[_spot_pix_num];
                _shifts = new int[_spot_pix_num];

                double spot_radius = (double)_exp_spot_dimension / 2;

                /** Haven't thought about this carefully for even-sized
                 *  spot models.
                 */
                double mid_point = (double)(model_diameter - 1) / 2;
                int mid_coord = model_diameter / 2;
                _model_half_dim = mid_coord;

                for (int i = 0; i < _spot_pix_num; i++)
                {
                        /* calculate distance of each pixel from the centre
                         * of the model spot
                         **/

                        int fast_coord = (i % model_diameter);
                        int slow_coord = (i / model_diameter);

                        double f_shift = fast_coord - mid_point;
                        double s_shift = slow_coord - mid_point;

                        double sqdist = pow(f_shift, 2) + pow(s_shift, 2);
                        double dist = sqrt(sqdist);

                        /* assume linear fall-off of spot height with
                         * increasing distance from centre
                         **/

                        double intensity = (spot_radius - dist) / spot_radius;
                        intensity *= BLOB_MAX_INTENSITY;

                        if (intensity < 0)
                        {
                                intensity = 0;
                        }

                        _ideal_spot[i] = intensity;

                        /** Shift corresponding to that on the complete
                         *  image. Need displacements in ints.
                         **/

                        int f_disp = fast_coord - mid_coord;
                        int s_disp = slow_coord - mid_coord;

                        _shifts[i] = s_disp * _pixels_per_row + f_disp;
                }

/* Delete this later */
#ifdef BLOB_DEBUG
                std::cout << "Ideal spot intensities (spot size"
                          << " " << _exp_spot_dimension << "):" << std::endl;

                for (int i = 0; i < model_diameter; i++)
                {
                        for (int j = 0; j < model_diameter; j++)
                        {
                                int index = i * model_diameter + j;

                                std::cout << _ideal_spot[index] << ", ";
                        }

                        std::cout << std::endl;
                }

                std::cout << std::endl;
                std::cout << "Spot shift positions:" << std::endl;

                for (int i = 0; i < model_diameter; i++)
                {
                        for (int j = 0; j < model_diameter; j++)
                        {
                                int index = i * model_diameter + j;

                                std::cout << _shifts[index] << ", ";
                        }

                        std::cout << std::endl;
                }

                std::cout << std::endl;

#endif // BLOB_DEBUG
        }

        template <typename T>
        int find_next_index(T *image, bool *mask, int start)
        {
                for (int i = start; i < _image_pix_num; i++)
                {
                        if (!mask[i])
                        {
                                continue;
                        }

                        if (image[i] > _global_threshold)
                        {
                                return i;
                        }
                }

                return -1;
        }

        template <typename T>
        int focus_on_maximum(T *image, bool *mask, int index)
        {
                T max = image[index];
                int best_idx = index;

                for (int i = 0; i < _spot_pix_num; i++)
                {
                        /* find index relative to current index
                         * to check the image value */
                        int test_idx = index + _shifts[i];

                        T value = image[test_idx];

                        if (value > max)
                        {
                                max = value;
                                best_idx = test_idx;
                        }
                }

                return best_idx;
        }

        bool is_viable_position(bool *mask, int index)
        {
                int slow_coord = (index / _pixels_per_row);
                int fast_coord = (index % _pixels_per_row);

                if ((slow_coord - _model_half_dim < 0) ||
                    (slow_coord + _model_half_dim >= _row_count))
                {
                        return false;
                }

                if ((fast_coord - _model_half_dim < 0) ||
                    (fast_coord + _model_half_dim >= _pixels_per_row))
                {
                        return false;
                }

                int masked_pixels = 0;
                int total_pixels = 0;

                for (int i = 0; i < _spot_pix_num; i++)
                {
                        /* find index relative to current index
                         * to check the image value */
                        int test_idx = index + _shifts[i];

                        masked_pixels += (mask[test_idx] == 0);
                        total_pixels++;
                }

                double ratio = (double)masked_pixels / (double)total_pixels;

                return (ratio < BLOB_MAX_MASKED_PIXELS);
        }

        template <typename T>
        double mean(T *values, int num)
        {
                T sum = 0;

                for (int i = 0; i < num; i++)
                {
                        sum += values[i];
                }

                sum /= num;

                return sum;
        }

        template <typename T>
        double correlation(T *xs, T *ys, int num)
        {
                double mean_x = mean(xs, num);
                double mean_y = mean(ys, num);

                double sum_x_y_rm_mean_x_y = 0;
                double sum_x_rm_mean_x_sq = 0;
                double sum_y_rm_mean_y_sq = 0;

                for (int i = 0; i < num; i++)
                {
                        T x = xs[i];
                        T y = ys[i];

                        sum_x_y_rm_mean_x_y += (x - mean_x) * (y - mean_y);
                        sum_x_rm_mean_x_sq += pow(x - mean_x, 2);
                        sum_y_rm_mean_y_sq += pow(y - mean_y, 2);
                }

                sum_x_y_rm_mean_x_y /= num;
                sum_x_rm_mean_x_sq /= num;
                sum_y_rm_mean_y_sq /= num;

                double r = sum_x_y_rm_mean_x_y;
                r /= sqrt(sum_x_rm_mean_x_sq * sum_y_rm_mean_y_sq);

                return r;
        }

        template <typename T>
        double blob_like_spot_score(T *image, bool *mask, int index)
        {
                /* We will use up to _spot_pix_num or less (indicated
                 * by count variable below */

                T *model_pix = new T[_spot_pix_num];
                T *data_pix = new T[_spot_pix_num];

                int count = 0;

                for (int i = 0; i < _spot_pix_num; i++)
                {
                        /* find index relative to current index
                         * to check the image value */
                        int test_idx = index + _shifts[i];

                        if (!mask[test_idx])
                        {
                                continue;
                        }

                        model_pix[count] = _ideal_spot[i];
                        data_pix[count] = image[test_idx];
                        count++;
                }

                double r = correlation(data_pix, model_pix, count);

                delete [] data_pix;

                return r;
        }

        void add_new_spot(bool *mask, bool *result, int index)
        {
                /* For now, really simple version. But in the future,
                 * more care to determine spot centroid. */

                /* Adding spot to 'result' */
                result[index] = 1;

                /* We do not want to consider surrounding pixels again */
                for (int i = 0; i < _spot_pix_num; i++)
                {
                        int test_idx = index + _shifts[i];
                        mask[test_idx] = 0;
                }
        }

        void clean_up()
        {
                delete [] _ideal_spot;
        }
public:
        BlobThresholdAlgorithm(int pixels_per_row, int row_count,
                               int exp_spot_dimension = 3,
                               double global_threshold = 100,
                               double min_blob_score = 0.7)
        {
                _pixels_per_row = pixels_per_row;
                _row_count = row_count;
                _image_pix_num = _pixels_per_row * _row_count;

                _min_blob_score = min_blob_score;
                _exp_spot_dimension = exp_spot_dimension;
                _global_threshold = global_threshold;

                calculate_ideal_spot();
        }

        template <typename T>
        void threshold(T *image, bool *mask, bool *result)
        {
                int index = 0;

                while (true)
                {
                        index = find_next_index(image, mask, index + 1);

                        if (index < 0)
                        {
                                break;
                        }

                        int max_index = index;

                        /* We would like to 'test' the spot at its centre */
                        max_index = focus_on_maximum(image, mask, index);

                        /* Is this spot on an edge of the image, or too
                         * many masked pixels in the spot window?
                         */
                        bool viable = is_viable_position(mask, index);

                        if (!viable)
                        {
                                continue;
                        }

                        double score;
                        score = blob_like_spot_score(image, mask, index);

                        if (score > _min_blob_score)
                        {
                                add_new_spot(mask, result, index);
                        }
                }
        }

        ~BlobThresholdAlgorithm()
        {
                clean_up();
        }
};

#endif // BLOBFINDER_H
