/*
 * blobfinder.h
 *
 *      An exploratory coding adventure into finding spot-shaped blobs on
 *      XFEL detectors.
 *
 *      Copyright (C) 2016 Helen Ginn
 *
 *      Author: Helen Ginn <helen@strubi.ax.ac.uk>
 *
 *      This code is distributed under the BSD license, a copy of which is
 *      included in the root directory of this package.
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

/* How much to pollute the initial 'ideal spot' with data */
#define IDEAL_SPOT_POLLUTION 0.5

/* Extra background pixels around ideal spot */
#define EXTRA_BACKGROUND_PIXELS 4

#ifdef BLOB_DEBUG
#include <iostream>
#endif

#include <cmath>

/* ------------------------------------------------------------------------
 * functions called under the core functions, still specialised (Level 3)
 * ------------------------------------------------------------------------*/

class BlobThresholdAlgorithm {
private:
  int _pixels_per_row;
  int _row_count;
  int _image_pix_num;
  int _exp_spot_dimension;
  int _model_half_dim;
  int _model_diameter;
  double _min_blob_score;
  int _num_passes;

  double *_ideal_spot;
  bool *_extraMask;
  int *_shifts;
  int _spot_pix_num;

  /* for now, hoping to remove this parameter */
  double _global_threshold;

  template <typename T>
  void update_ideal_spot(T *image, bool *mask, bool *result) {
    double *sum_found = new double[_spot_pix_num];
    int spot_count = 0;

    for (int i = 0; i < _image_pix_num; i++) {
      bool is_spot = result[i];

      if (!is_spot) {
        continue;
      }

      for (int j = 0; j < _spot_pix_num; j++) {
        int test_idx = i + _shifts[j];

        T value = image[test_idx];
        sum_found[j] += value;
      }

      spot_count++;
    }

    double ideal_spot_sum = 0;
    double pollute_spot_sum = 0;

    /* Find the total sum of both current ideal spot and the
     * new 'polluting' data from our image */

    for (int i = 0; i < _spot_pix_num; i++) {
      ideal_spot_sum += _ideal_spot[i];
      pollute_spot_sum += sum_found[i];
    }

    double normalise_factor = ideal_spot_sum / pollute_spot_sum;

    double pollution = IDEAL_SPOT_POLLUTION;
    double clean = 1 - pollution;

    for (int i = 0; i < _spot_pix_num; i++) {
      sum_found[i] *= normalise_factor;

      _ideal_spot[i] *= clean;
      _ideal_spot[i] += sum_found[i] * pollution;
    }

    delete[] sum_found;

    //      std::cout << "Updated spot values: " << std::endl;
    //      print_current_spot();
  }

  void set_ideal_spot(double *new_spot, int pix_num) {
    if (pix_num != _spot_pix_num) {
      std::cout << "Warning: ideal spot pixels"
                << " not same number." << std::endl;
      std::cout << "Not continuing." << std::endl;

      return;
    }

    memcpy(_ideal_spot, new_spot, pix_num * sizeof(double));
  }

  /**
   * Calculates what an 'ideal' spot looks like and stores this in
   * the member _ideal_spot[].
   * Also calculates the 'shifts' from 0 assuming the centre of an
   * odd-pixel diameter spot is in the centre, and for an even-diameter
   * spot, displaced to the bottom right slightly.
   **/
  void calculate_ideal_spot() {
    /* add two "background" pixels for comparison */
    _model_diameter = _exp_spot_dimension + EXTRA_BACKGROUND_PIXELS;

    _spot_pix_num = pow((double)_model_diameter, 2);

    _ideal_spot = new double[_spot_pix_num];
    _shifts = new int[_spot_pix_num];

    double spot_radius = (double)_exp_spot_dimension / 2;

    /** Haven't thought about this carefully for even-sized
     *      spot models.
     */
    double mid_point = (double)(_model_diameter - 1) / 2;
    int mid_coord = _model_diameter / 2;
    _model_half_dim = mid_coord;

    for (int i = 0; i < _spot_pix_num; i++) {
      /* calculate distance of each pixel from the centre
       * of the model spot
       **/

      int fast_coord = (i % _model_diameter);
      int slow_coord = (i / _model_diameter);

      double f_shift = fast_coord - mid_point;
      double s_shift = slow_coord - mid_point;

      double sqdist = pow(f_shift, 2) + pow(s_shift, 2);
      double dist = sqrt(sqdist);

      /* assume linear fall-off of spot height with
       * increasing distance from centre
       **/

      double intensity = (spot_radius - dist) / spot_radius;
      intensity *= BLOB_MAX_INTENSITY;

      if (intensity < 0) {
        intensity = 0;
      }

      _ideal_spot[i] = intensity;

      /** Shift corresponding to that on the complete
       *      image. Need displacements in ints.
       **/

      int f_disp = fast_coord - mid_coord;
      int s_disp = slow_coord - mid_coord;

      _shifts[i] = s_disp * _pixels_per_row + f_disp;
    }

    print_current_spot();
  }

  /* Delete this later */

  void print_current_spot() {
#ifdef BLOB_DEBUG
    std::cout << "Ideal spot intensities (spot size"
              << " " << _exp_spot_dimension << "):" << std::endl;

    for (int i = 0; i < _model_diameter; i++) {
      for (int j = 0; j < _model_diameter; j++) {
        int index = i * _model_diameter + j;

        std::cout << _ideal_spot[index] << ", ";
      }

      std::cout << std::endl;
    }
#endif  // BLOB_DEBUG
  }

  bool isMasked(bool *mask, int i) {
    bool isNotMasked = mask[i] * _extraMask[i];

    return !isNotMasked;
  }

  template <typename T>
  int find_next_index(T *image, bool *mask, int start) {
    for (int i = start; i < _image_pix_num; i++) {
      if (!mask[i]) {
        continue;
      }

      if (image[i] > _global_threshold) {
        return i;
      }
    }

    return -1;
  }

  template <typename T>
  int focus_on_maximum(T *image, bool *mask, int index) {
    T max = image[index];
    int best_idx = index;

    for (int i = 0; i < _spot_pix_num; i++) {
      /* find index relative to current index
       * to check the image value */
      int test_idx = index + _shifts[i];

      T value = image[test_idx];

      if (value > max) {
        max = value;
        best_idx = test_idx;
      }
    }

    return best_idx;
  }

  bool is_viable_position(bool *mask, int index) {
    int slow_coord = (index / _pixels_per_row);
    int fast_coord = (index % _pixels_per_row);

    if ((slow_coord - _model_half_dim < 0)
        || (slow_coord + _model_half_dim >= _row_count)) {
      return false;
    }

    if ((fast_coord - _model_half_dim < 0)
        || (fast_coord + _model_half_dim >= _pixels_per_row)) {
      return false;
    }

    int masked_pixels = 0;
    int total_pixels = 0;

    for (int i = 0; i < _spot_pix_num; i++) {
      /* find index relative to current index
       * to check the image value */
      int test_idx = index + _shifts[i];
      bool masked = isMasked(mask, test_idx);

      masked_pixels += (masked);
      total_pixels++;
    }

    double ratio = (double)masked_pixels / (double)total_pixels;

    return (ratio < BLOB_MAX_MASKED_PIXELS);
  }

  template <typename T>
  double mean(T *values, int num) {
    T sum = 0;

    for (int i = 0; i < num; i++) {
      sum += values[i];
    }

    sum /= num;

    return sum;
  }

  template <typename T>
  double correlation(T *xs, T *ys, int num) {
    double mean_x = mean(xs, num);
    double mean_y = mean(ys, num);

    double sum_x_y_rm_mean_x_y = 0;
    double sum_x_rm_mean_x_sq = 0;
    double sum_y_rm_mean_y_sq = 0;

    for (int i = 0; i < num; i++) {
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
  double blob_like_spot_score(T *image, bool *mask, int index) {
    /* Is this spot on an edge of the image, or too
     * many masked pixels in the spot window?
     */
    bool viable = is_viable_position(mask, index);

    if (!viable) {
      return 0;
    }

    /* We will use up to _spot_pix_num or less (indicated
     * by count variable below */

    T *model_pix = new T[_spot_pix_num];
    T *data_pix = new T[_spot_pix_num];

    int count = 0;

    for (int i = 0; i < _spot_pix_num; i++) {
      /* find index relative to current index
       * to check the image value */
      int test_idx = index + _shifts[i];

      if (!mask[test_idx]) {
        continue;
      }

      model_pix[count] = _ideal_spot[i];
      data_pix[count] = image[test_idx];
      count++;
    }

    double r = correlation(data_pix, model_pix, count);

    delete[] data_pix;

    return r;
  }

  void add_new_spot(bool *mask, bool *result, int index) {
    /* For now, really simple version. But in the future,
     * more care to determine spot centroid. */

    /* Adding spot to 'result' */
    result[index] = 1;

    /* We do not want to consider surrounding pixels again */
    for (int i = 0; i < _spot_pix_num; i++) {
      int test_idx = index + _shifts[i];
      _extraMask[test_idx] = 0;
    }
  }

  void clean_up() {
    delete[] _extraMask;
    delete[] _ideal_spot;
  }

public:
  BlobThresholdAlgorithm(int pixels_per_row,
                         int row_count,
                         int exp_spot_dimension = 3,
                         double global_threshold = 100,
                         double min_blob_score = 0.7,
                         int num_passes = 0) {
    _pixels_per_row = pixels_per_row;
    _row_count = row_count;
    _image_pix_num = _pixels_per_row * _row_count;

    _min_blob_score = min_blob_score;
    _exp_spot_dimension = exp_spot_dimension;
    _global_threshold = global_threshold;
    _num_passes = num_passes;

    calculate_ideal_spot();
  }

  template <typename T>
  void correlation_image(T *image, bool *mask, double *result) {
    int index = 0;

    while (true) {
      index = find_next_index(image, mask, index + 1);

      if (index < 0) {
        break;
      }

      int max_index = index;

      /* We would like to 'test' the spot at its centre */
      max_index = focus_on_maximum(image, mask, index);

      double score;
      score = blob_like_spot_score(image, mask, index);
      result[index] = score;
    }
  }

  template <typename T>
  void find_results(T *image, bool *mask, bool *result) {
    int index = 0;
    int spot_count = 0;

    while (true) {
      index = find_next_index(image, mask, index + 1);

      if (index < 0) {
        break;
      }

      int max_index = index;

      /* We would like to 'test' the spot at its centre */
      max_index = focus_on_maximum(image, mask, index);

      double score;
      score = blob_like_spot_score(image, mask, index);

      if (score > _min_blob_score) {
        spot_count++;
        add_new_spot(mask, result, index);
      }
    }

    std::cout << "Spot count after pass: " << spot_count << std::endl;
  }

  void reset_results(bool *result) {
    memset(result, 0, sizeof(bool) * _image_pix_num);
    memset(_extraMask, 1, sizeof(bool) * _image_pix_num);
  }

  template <typename T>
  void threshold(T *image, bool *mask, bool *result) {
    _extraMask = new bool[_image_pix_num];
    reset_results(result);

    find_results(image, mask, result);

    for (int i = 0; i < _num_passes; i++) {
      update_ideal_spot(image, mask, result);
      reset_results(result);
      find_results(image, mask, result);
    }
  }

  ~BlobThresholdAlgorithm() {
    clean_up();
  }
};

#endif  // BLOBFINDER_H
