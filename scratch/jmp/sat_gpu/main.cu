#include <cstdio>
#include <iostream>
#include <cassert>
#include <cudpp.h>
#include <ctime>
#include <sys/time.h>

#define GPUERROR(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) {
      exit(code);
    }
  }
}


template <typename T>
__global__ void transpose_kernel(T *dst, const T *src, std::size_t width, std::size_t height) {

  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  std::size_t j = blockIdx.y * blockDim.y + threadIdx.y;
  std::size_t k1 = i + j * width;
  std::size_t k2 = j + i * height;
  if (i < width && j < height) {
    dst[k2] = src[k1];
  }
}


template <typename T>
CUDPPConfiguration sat_scan_configuration();

template <>
CUDPPConfiguration sat_scan_configuration<int>() {
  CUDPPConfiguration config = { 
    CUDPP_SCAN, 
    CUDPP_ADD, 
    CUDPP_INT,
    CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE 
  };
  return config;
}

template <>
CUDPPConfiguration sat_scan_configuration<float>() {
  CUDPPConfiguration config = { 
    CUDPP_SCAN, 
    CUDPP_ADD, 
    CUDPP_FLOAT,
    CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE 
  };
  return config;
}

template <>
CUDPPConfiguration sat_scan_configuration<double>() {
  CUDPPConfiguration config = { 
    CUDPP_SCAN, 
    CUDPP_ADD, 
    CUDPP_DOUBLE,
    CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE 
  };
  return config;
}

template <typename T>
class SummedAreaTableGpu {
public: 

  SummedAreaTableGpu(std::size_t width, std::size_t height) {
   
    // Save width and height
    width_ = width;
    height_ = height;

    // Initialise the cudpp handle
    cudppCreate(&cudpp_handle_);

    // Setup the scan configuration
    CUDPPConfiguration config = sat_scan_configuration<T>();

    // Create the scan plans
    cudppPlan(cudpp_handle_, &scan_rows_, config, width, height, width);
    cudppPlan(cudpp_handle_, &scan_cols_, config, height, width, height);

    // Allocate the buffer
    GPUERROR(cudaMalloc((void **)&buffer_, width*height*sizeof(T)));
  }

  ~SummedAreaTableGpu() {
   
    // Free the buffer
    cudaFree(buffer_);

    // Free the scan plans
    cudppDestroyPlan(scan_rows_);
    cudppDestroyPlan(scan_cols_);

    // Free the cudpp handle
    cudppDestroy(cudpp_handle_);
  }

  void operator()(T *dst, const T *src) {

    // Scan the rows of the image
    cudppMultiScan(scan_rows_, buffer_, src, width_, height_);

    // Transpose the image
    transpose(dst, buffer_, width_, height_);

    // Scan the columns of the image
    cudppMultiScan(scan_cols_, buffer_, dst, height_, width_);

    // Transpose the image
    transpose(dst, buffer_, height_, width_);
  }

private:

  void transpose(T *dst, const T *src, std::size_t width, std::size_t height) {
    
    // Set the block and grid sizes
    dim3 block(16, 16, 1);
    dim3 grid(
      width % block.x == 0 ? width / block.x : width / block.x + 1, 
      height % block.y == 0 ? height / block.y : height / block.y + 1,
      1);

    // Transpose the image
    transpose_kernel<T><<<grid, block>>>(dst, src, width, height);
    GPUERROR(cudaPeekAtLastError());
  }

  std::size_t width_;
  std::size_t height_;
  CUDPPHandle cudpp_handle_;
  CUDPPHandle scan_rows_;
  CUDPPHandle scan_cols_;
  T *buffer_;
};


template <typename T>
void sat_cpu(T *dst, const T *src, std::size_t width, std::size_t height) {
  for (std::size_t j = 0; j < height; ++j) {
    for (std::size_t i = 0; i < width; ++i) {
      T I10 = (j > 0 ? dst[i+(j-1)*width] : 0);
      T I01 = (i > 0 ? dst[(i-1)+j*width] : 0);
      T I11 = (i > 0 && j > 0 ? dst[(i-1)+(j-1)*width] : 0);
      dst[i+j*width] = src[i+j*width] + I10 + I01 - I11;
    }
  }
}

template <typename T>
void check_sat_results(const T *res1, const T *res2, std::size_t width, std::size_t height) {
  std::cout << " Checking results are correct" << std::endl;
  double EPS = 1e-3;
  bool success = true;
  for (std::size_t i = 0; i < width * height; ++i) {
    double d = (double)res1[i] - (double)res2[i];
    if (std::abs(d / (double)res1[i]) > EPS) {
      success = false;
      break;
    }
  }
  std::cout << "  Test " << (success ? "PASSED" : "FAILED") << std::endl;
  assert(success == true);
}

typedef unsigned long long timestamp_t;

inline
timestamp_t get_timestamp() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_usec + (timestamp_t)t.tv_sec * 1000000;
}

template <typename T>
timestamp_t test_sat_cpu(T *dst, const T *src, std::size_t width, 
    std::size_t height, std::size_t niter) {
  std::cout << " Testing CPU SAT" << std::endl;
  timestamp_t st = get_timestamp();
  for (std::size_t i = 0; i < niter; ++i) {
    sat_cpu<T>(dst, src, width, height);
  }
  timestamp_t ft = get_timestamp();
  std::cout << "  Time taken: " << (ft - st) / 1000000.0 << " seconds" << std::endl;
  return ft - st;
}

template <typename T>
timestamp_t test_sat_gpu(T *dst, T *src, std::size_t width, 
    std::size_t height, std::size_t niter) {
  std::cout << " Testing GPU SAT" << std::endl;
  SummedAreaTableGpu<T> summed_area_table(width, height);
  timestamp_t st = get_timestamp();
  for (std::size_t i = 0; i < niter; ++i) {
    summed_area_table(dst, src);
  }
  cudaDeviceSynchronize();
  timestamp_t ft = get_timestamp();
  std::cout << "  Time taken: " << (ft - st) / 1000000.0 << " seconds" << std::endl;
  return ft - st;
}


template <typename T>
void test_sat(std::size_t width, std::size_t height, std::size_t niter) {

  // Declare the arrays we need
  T *gpu_src = 0, *gpu_dst = 0;
  T *cpu_src = new T[width * height];
  T *cpu_dst = new T[width * height];
  T *cpu_dst2 = new T[width * height];

  // Init the image to process
  for (std::size_t i = 0; i < width * height; ++i) {
    cpu_src[i] = rand();
  }

  // Allocate the gpu arrays and copy the data
  GPUERROR(cudaMalloc((void **)&gpu_src, width*height*sizeof(T)));
  GPUERROR(cudaMalloc((void **)&gpu_dst, width*height*sizeof(T)));
  GPUERROR(cudaMemcpy((void *)gpu_src, (void *)cpu_src, 
    width*height*sizeof(T),  cudaMemcpyHostToDevice));

  // Test the summed area table
  timestamp_t t1 = test_sat_cpu<T>(cpu_dst, cpu_src, width, height, niter);
  timestamp_t t2 = test_sat_gpu<T>(gpu_dst, gpu_src, width, height, niter);
  std::cout << "GPU Speedup = " << (double)t1 / (double) t2 << " X" << std::endl;

  // Check the results
  GPUERROR(cudaMemcpy((void *)cpu_dst2, (void *)gpu_dst, width*height*sizeof(T), cudaMemcpyDeviceToHost));
  check_sat_results(cpu_dst, cpu_dst2, width, height);

  // Delete the arrays 
  cudaFree(gpu_src);
  cudaFree(gpu_dst);
  delete[] cpu_src;
  delete[] cpu_dst;
  delete[] cpu_dst2;
}

int main() {

  test_sat<int>(1024, 1024, 100);
  test_sat<float>(1024, 1024, 100);
  /*test_sat<double>(1024, 1024, 100);*/


  return 0;
}
