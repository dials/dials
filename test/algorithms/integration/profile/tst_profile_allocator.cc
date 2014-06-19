
#include <dials/algorithms/integration/profile/profile_allocator.h>

using namespace dials;
using namespace dials::algorithms;


int main() {

  // Setup the sizes
  std::size_t nrefl = 100000;
  std::size_t size = 10;
  std::size_t num = nrefl * size;

  // Allocate arrays
  af::shared<int> frame(num);
  af::shared<int4> bbox(num);
  std::size_t radius = 10;

  // Fill arrays
  for (std::size_t i = 0, k = 0; i < nrefl; ++i) {
    int z0 = rand() % 1000;
    int x0 = rand() % 1000;
    int y0 = rand() % 1000;
    int nx = 3 + rand() % 10;
    int ny = 3 + rand() % 10;
    int x1 = x0 + nx;
    int y1 = y0 + ny;
    for (std::size_t j = 0; j < size; ++j, ++k) {
      frame[k] = z0 + j;
      bbox[k][0] = x0;
      bbox[k][1] = x1;
      bbox[k][2] = y0;
      bbox[k][3] = y1;
    }
  }

  // Create the profile allocator
  ProfileAllocator allocator(frame.const_ref(), bbox.const_ref(), radius);

  // Check the maximum size if correct
  std::size_t max_size = 0;
  for (std::size_t i = 0; i < bbox.size(); ++i) {
    std::size_t size = (bbox[i][1]-bbox[i][0])*(bbox[i][3]-bbox[i][2]);
    max_size = std::max(max_size, size);
  }
  DIALS_ASSERT(max_size == allocator.max_size());
  std::cout << "OK" << std::endl;

  // Set the number of spots on each frame
  std::vector<std::size_t> count(1000+10, 0);
  for (std::size_t i = 0; i < frame.size(); ++i) {
    count[frame[i]]++;
  }

  // Check the maximum number is ok
  std::size_t max_image = 0;
  for (std::size_t i = 0; i < count.size(); ++i) {
    max_image = std::max(max_image, count[i]);
  }
  DIALS_ASSERT(max_image == allocator.max_image());
  std::cout << "OK" << std::endl;

  // Check the maxmimum number is ok
  std::size_t max_num = 0;
  for (int i = 0; i < count.size(); ++i) {
    std::size_t num = 0;
    for (int j = i - radius; j <= i+radius; ++j) {
      if (j >= 0 && j < count.size()) {
        num += count[j];
      }
    }
    max_num = std::max(max_num, num);
  }
  DIALS_ASSERT(max_num == allocator.max_num());
  std::cout << "OK" << std::endl;


  // Check that getting the profile is ok
  std::vector<std::size_t> index(frame.size());
  for (std::size_t i = 0; i < index.size(); ++i) {
    index[i] = i;
  }
  std::sort(index.begin(), index.end(), ProfileAllocator::sort_by_z(frame.const_ref()));
  std::size_t mini = 0, maxi = 0;
  for (int f = 0; f < 1000+10; ++f) {
    while (mini < index.size() && frame[index[mini]] < (f - (int)(2*radius-1))) {
      int *data = allocator.data(index[mini]);
      for (std::size_t i = 0; i < max_size; ++i) {
        DIALS_ASSERT(data[i] == index[mini]);
      }
      allocator.free(index[mini]);
      mini++;
    }
    while (maxi < index.size() && frame[index[maxi]] == f) {
      allocator.hold(index[maxi]);
      int * data = allocator.data(index[maxi]);
      for (std::size_t i = 0; i < max_size; ++i) {
        data[i] = index[maxi];
      }
      maxi++;
    }
  }
  std::cout << "OK" << std::endl;

  return 0;
}
