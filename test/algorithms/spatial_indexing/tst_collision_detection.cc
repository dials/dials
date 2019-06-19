#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <deque>
#include <dials/algorithms/spatial_indexing/detect_collisions.h>

struct Box {
  int x0, y0, x1, y1;
  Box() {}
  Box(int x0_, int y0_, int x1_, int y1_) : x0(x0_), y0(y0_), x1(x1_), y1(y1_) {}
};

struct Box3d {
  int x0, y0, z0, x1, y1, z1;
  Box3d() {}
  Box3d(int x0_, int y0_, int z0_, int x1_, int y1_, int z1_)
      : x0(x0_), y0(y0_), z0(z0_), x1(x1_), y1(y1_), z1(z1_) {}
};

Box random_box(Box bounds, int min_size, int max_size) {
  int x = bounds.x0 + rand() % (bounds.x1 - bounds.x0);
  int y = bounds.y0 + rand() % (bounds.y1 - bounds.y0);
  int w = min_size + rand() % (max_size - min_size);
  int h = min_size + rand() % (max_size - min_size);
  return Box(x, y, x + w, y + h);
}

Box3d random_box3d(Box3d bounds, int min_size, int max_size) {
  int x = bounds.x0 + rand() % (bounds.x1 - bounds.x0);
  int y = bounds.y0 + rand() % (bounds.y1 - bounds.y0);
  int z = bounds.z0 + rand() % (bounds.z1 - bounds.z0);
  int w = min_size + rand() % (max_size - min_size);
  int h = min_size + rand() % (max_size - min_size);
  int d = min_size + rand() % (max_size - min_size);
  return Box3d(x, y, z, x + w, y + h, z + d);
}

// std::ostream& operator<<(std::ostream &os, const Box &box) {
//  os << "(" << box.x0 << ", " << box.x1 << ", " << box.y0 << ", " << box.y1 << ")";
//  return os;
//}

// std::string range(int minx, int maxx, int miny, int maxy) {
//  std::stringstream ss;
//  ss << "(" << minx << ", " << maxx << ", " << miny << ", " << maxy << ")";
//  return ss.str();
//}

using dials::algorithms::detect_collisions2d;
using dials::algorithms::detect_collisions3d;
using dials::algorithms::get_maximum_bound;
using dials::algorithms::get_minimum_bound;

namespace dials { namespace algorithms {

  // Helper functions needed for 2D collision detection
  template <>
  double get_minimum_bound<0, Box>(const Box &b) {
    return b.x0;
  }
  template <>
  double get_minimum_bound<1, Box>(const Box &b) {
    return b.y0;
  }
  template <>
  double get_maximum_bound<0, Box>(const Box &b) {
    return b.x1;
  }
  template <>
  double get_maximum_bound<1, Box>(const Box &b) {
    return b.y1;
  }

  // Helper functions needed for 3D collision detection
  template <>
  double get_minimum_bound<0, Box3d>(const Box3d &b) {
    return b.x0;
  }
  template <>
  double get_minimum_bound<1, Box3d>(const Box3d &b) {
    return b.y0;
  }
  template <>
  double get_minimum_bound<2, Box3d>(const Box3d &b) {
    return b.z0;
  }
  template <>
  double get_maximum_bound<0, Box3d>(const Box3d &b) {
    return b.x1;
  }
  template <>
  double get_maximum_bound<1, Box3d>(const Box3d &b) {
    return b.y1;
  }
  template <>
  double get_maximum_bound<2, Box3d>(const Box3d &b) {
    return b.z1;
  }

}}  // namespace dials::algorithms

void tst_detect_2d() {
  int num = 10000;
  std::vector<Box> data(num);
  std::deque<std::pair<int, int> > collisions1;
  std::deque<std::pair<int, int> > collisions2;

  Box bounds(0, 0, 20000, 20000);

  // Create a number of random boxes
  for (std::size_t i = 0; i < num; ++i) {
    data[i] = random_box(bounds, 3, 8);
  }

  // Do the collision check
  //  clock_t st = clock();
  detect_collisions2d(data.begin(), data.end(), collisions1);
  //  std::cout << ((float)(clock() - st)) / CLOCKS_PER_SEC << "\n";
  //  std::cout << collisions1.size() << "\n";
  //
  // Do a brute force check to see if we get the correct results.
  for (std::size_t j = 0; j < num - 1; ++j) {
    int jx0 = data[j].x0;
    int jx1 = data[j].x1;
    int jy0 = data[j].y0;
    int jy1 = data[j].y1;
    for (std::size_t i = j + 1; i < num; ++i) {
      int ix0 = data[i].x0;
      int ix1 = data[i].x1;
      int iy0 = data[i].y0;
      int iy1 = data[i].y1;
      if (!(ix0 >= jx1 || jx0 >= ix1 || iy0 >= jy1 || jy0 >= iy1)) {
        collisions2.push_back(std::pair<int, int>(i, j));
      }
    }
  }

  // Check the sizes are equal
  assert(collisions2.size() == collisions1.size());

  // Test passed
  std::cout << "OK" << std::endl;
}

void tst_detect_3d() {
  int num = 10000;
  std::vector<Box3d> data(num);
  std::deque<std::pair<int, int> > collisions1;
  std::deque<std::pair<int, int> > collisions2;

  Box3d bounds(0, 0, 0, 512, 512, 512);

  // Create a load of random boxes
  for (std::size_t i = 0; i < num; ++i) {
    data[i] = random_box3d(bounds, 3, 8);
  }

  // Do the collision check
  detect_collisions3d(data.begin(), data.end(), collisions1);

  // Do a brute force check
  for (std::size_t j = 0; j < num - 1; ++j) {
    int jx0 = data[j].x0;
    int jx1 = data[j].x1;
    int jy0 = data[j].y0;
    int jy1 = data[j].y1;
    int jz0 = data[j].z0;
    int jz1 = data[j].z1;
    for (std::size_t i = j + 1; i < num; ++i) {
      int ix0 = data[i].x0;
      int ix1 = data[i].x1;
      int iy0 = data[i].y0;
      int iy1 = data[i].y1;
      int iz0 = data[i].z0;
      int iz1 = data[i].z1;

      if (!(ix0 >= jx1 || jx0 >= ix1 || iy0 >= jy1 || jy0 >= iy1 || iz0 >= jz1
            || jz0 >= iz1)) {
        collisions2.push_back(std::pair<int, int>(i, j));
      }
    }
  }

  // Check the sizes are equal
  assert(collisions2.size() == collisions1.size());

  // Test passed
  std::cout << "OK" << std::endl;
}

int main(int argc, char const *argv[]) {
  tst_detect_2d();
  tst_detect_3d();

  return 0;
}
