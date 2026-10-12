/*
 * bench_connected_components.cc
 *
 * Compare connected component labelling of pixel lists before and after the
 * switch from boost::adjacency_list + boost::connected_components to the
 * DisjointSets union-find (#3293).
 *
 * The "Before" benchmarks run a verbatim copy of the old implementation,
 * kept below in the `before` namespace. The "After" benchmarks call the
 * current DIALS headers. Both are given identical input, and the labels they
 * produce are checked to be identical before any timing is done.
 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef __GLIBC__
#include <malloc.h>
#endif

#include <benchmark/benchmark.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/unordered_map.hpp>

#include <dials/algorithms/image/connected_components/connected_components.h>
#include <dials/model/data/pixel_list.h>

using dials::algorithms::LabelPixels;
using dials::model::PixelList;
using dials::model::PixelListLabeller;
using scitbx::vec3;
using scitbx::af::int2;
using scitbx::af::int3;
namespace af = scitbx::af;

// ---------------------------------------------------------------------------
// The previous implementation, copied from before the union-find change
// ---------------------------------------------------------------------------
namespace before {

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
  AdjacencyList;

/** PixelListLabeller::labels_3d and labels_2d, before */
af::shared<int> pixel_list_labels(const af::const_ref<vec3<int> >& coords_,
                                  int2 size_,
                                  int last_frame_,
                                  bool threed) {
  if (coords_.size() == 0) {
    return af::shared<int>();
  }

  // Calculate the coordinate indices
  for (std::size_t i = 0; i < coords_.size(); ++i) {
    if (i > 0) {
      DIALS_ASSERT(dials::model::detail::lessthan(coords_[i - 1], coords_[i]));
    }
  }

  // Create a graph of coordinates
  AdjacencyList graph(coords_.size());
  std::size_t i1 = 0, i2 = 0, i3 = 0;
  for (; i1 < coords_.size() - 1; ++i1) {
    vec3<int> a0 = coords_[i1];
    vec3<int> a1(a0[0], a0[1], a0[2] + 1);
    vec3<int> a2(a0[0], a0[1] + 1, a0[2]);
    vec3<int> a3(a0[0] + 1, a0[1], a0[2]);
    if (coords_[i1 + 1] == a1) {
      boost::add_edge(i1, i1 + 1, graph);
    }
    if (a0[1] < size_[0] - 1) {
      for (; dials::model::detail::lessthan(coords_[i2], a2) && i2 < coords_.size() - 1;
           ++i2)
        ;
      if (coords_[i2] == a2) {
        boost::add_edge(i1, i2, graph);
      }
    }
    if (threed && a0[0] < last_frame_ - 1) {
      if (i2 > i3) i3 = i2;
      for (; dials::model::detail::lessthan(coords_[i3], a3) && i3 < coords_.size() - 1;
           ++i3)
        ;
      if (coords_[i3] == a3) {
        boost::add_edge(i1, i3, graph);
      }
    }
  }

  // Do the connected components
  af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
  DIALS_ASSERT(labels.size() == coords_.size());
  int num = boost::connected_components(graph, &labels[0]);
  DIALS_ASSERT(num <= labels.size());
  return labels;
}

struct Vec3IntHash {
  vec3<int> size_;
  Vec3IntHash(vec3<int> size) : size_(size) {}
  std::size_t operator()(vec3<int> const& x) const {
    boost::hash<int> hasher;
    return hasher(x[0] + x[1] * size_[2] + x[2] * size_[1] * size_[2]);
  }
};

/** LabelPixels::labels, before. coords_ are (z, y, x) */
af::shared<int> label_pixels_labels(const af::const_ref<vec3<int> >& coords_,
                                    int3 size_) {
  // Create a hash table of the points
  AdjacencyList graph(coords_.size());
  boost::unordered_map<vec3<int>, int, Vec3IntHash> grid(coords_.size(),
                                                         Vec3IntHash(size_));
  for (std::size_t i = 0; i < coords_.size(); ++i) {
    DIALS_ASSERT(grid[coords_[i]] == 0);
    grid[coords_[i]] = i + 1;
  }

  // For each point check the pixels to the left in all three dimensions
  // and if they are in the list of pixels then add the edges.
  for (std::size_t i = 0; i < coords_.size(); ++i) {
    int j;
    j = grid[vec3<int>(coords_[i][0] - 1, coords_[i][1], coords_[i][2])];
    if (j != 0) boost::add_edge(i, j - 1, graph);
    j = grid[vec3<int>(coords_[i][0], coords_[i][1] - 1, coords_[i][2])];
    if (j != 0) boost::add_edge(i, j - 1, graph);
    j = grid[vec3<int>(coords_[i][0], coords_[i][1], coords_[i][2] - 1)];
    if (j != 0) boost::add_edge(i, j - 1, graph);
  }

  // Do the connected components
  af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
  int num = boost::connected_components(graph, &labels[0]);
  DIALS_ASSERT(num <= labels.size());
  return labels;
}

}  // namespace before

// ---------------------------------------------------------------------------
// Contrived inputs. Each is a stack of `frames` square images of side `side`.
// ---------------------------------------------------------------------------
namespace {

enum Pattern {
  /** ~2% of pixels set at random: mostly isolated pixels, few edges */
  Sparse,
  /** 3x3x3 spots on a grid, with gaps between them: many small components */
  Spots,
  /** Every pixel set: one component with the maximum number of edges */
  Solid,
  /**
   * A single-pixel-wide boustrophedon path through each frame, joined
   * between frames: one long, thin component (deep union-find trees)
   */
  Serpentine,
};

const char* pattern_name(Pattern p) {
  switch (p) {
  case Sparse:
    return "Sparse";
  case Spots:
    return "Spots";
  case Solid:
    return "Solid";
  case Serpentine:
    return "Serpentine";
  }
  return "?";
}

bool is_set(Pattern p, int k, int j, int i, int side, std::mt19937& rng) {
  switch (p) {
  case Sparse:
    // 2% of the full range of mt19937
    return rng() < 85899346u;
  case Spots:
    // 3 pixels on, 5 off, in every dimension
    return (k % 8 < 3) && (j % 8 < 3) && (i % 8 < 3);
  case Solid:
    return true;
  case Serpentine:
    // Even rows are full; odd rows have one pixel at alternating ends
    if (j % 2 == 0) return true;
    return ((j / 2) % 2 == 0) ? (i == side - 1) : (i == 0);
  }
  return false;
}

/** The mask for frame k of the pattern */
af::versa<bool, af::c_grid<2> > make_mask(Pattern p,
                                          int side,
                                          int k,
                                          std::mt19937& rng) {
  af::versa<bool, af::c_grid<2> > mask(af::c_grid<2>(side, side), false);
  for (int j = 0; j < side; ++j) {
    for (int i = 0; i < side; ++i) {
      mask(j, i) = is_set(p, k, j, i, side, rng);
    }
  }
  return mask;
}

PixelListLabeller make_labeller(Pattern p, int side, int frames) {
  std::mt19937 rng(42);
  af::versa<int, af::c_grid<2> > image(af::c_grid<2>(side, side), 1);
  PixelListLabeller labeller;
  for (int k = 0; k < frames; ++k) {
    af::versa<bool, af::c_grid<2> > mask = make_mask(p, side, k, rng);
    labeller.add(PixelList(k, image.const_ref(), mask.const_ref()));
  }
  return labeller;
}

LabelPixels make_label_pixels(Pattern p, int side, int frames) {
  std::mt19937 rng(42);
  // LabelPixels takes (x, y, z) coordinates. Give them in a scrambled order
  // since, unlike the pixel list, it does not rely on them being sorted.
  af::shared<vec3<int> > xyz;
  for (int k = 0; k < frames; ++k) {
    af::versa<bool, af::c_grid<2> > mask = make_mask(p, side, k, rng);
    for (int j = 0; j < side; ++j) {
      for (int i = 0; i < side; ++i) {
        if (mask(j, i)) xyz.push_back(vec3<int>(i, j, k));
      }
    }
  }
  std::shuffle(xyz.begin(), xyz.end(), std::mt19937(7));
  af::shared<int> values(xyz.size(), 1);
  LabelPixels label_pixels(int3(frames, side, side));
  label_pixels.add_pixels(values.const_ref(), xyz.const_ref());
  return label_pixels;
}

void check_same(const af::shared<int>& a, const af::shared<int>& b) {
  if (a.size() != b.size() || !std::equal(a.begin(), a.end(), b.begin())) {
    std::cerr << "Labels from before and after do not match" << std::endl;
    std::abort();
  }
}

/**
 * Hand freed memory back to the system. The Before labelling frees millions
 * of small graph nodes, which leaves glibc's heap in a state that makes
 * later large allocations several times slower. Without this, the After
 * benchmarks would be charged for the check against Before.
 */
void release_heap() {
#ifdef __GLIBC__
  malloc_trim(0);
#endif
}

void set_counters(benchmark::State& state,
                  std::size_t num_pixels,
                  const af::shared<int>& labels) {
  int num_labels = 0;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    num_labels = std::max(num_labels, labels[i] + 1);
  }
  state.SetItemsProcessed(state.iterations() * num_pixels);
  state.counters["pixels"] = num_pixels;
  state.counters["spots"] = num_labels;
}

// -------------------------------------------------------------------------
// PixelListLabeller (dials.find_spots)
// -------------------------------------------------------------------------

void BM_PixelList(benchmark::State& state, Pattern p, bool threed, bool after) {
  const int side = state.range(0);
  const int frames = state.range(1);
  PixelListLabeller labeller = make_labeller(p, side, frames);
  af::shared<vec3<int> > coords = labeller.coords();

  af::shared<int> expected = threed ? labeller.labels_3d() : labeller.labels_2d();
  check_same(expected,
             before::pixel_list_labels(
               coords.const_ref(), labeller.size(), labeller.last_frame(), threed));
  release_heap();

  auto run = [&]() {
    return after
             ? (threed ? labeller.labels_3d() : labeller.labels_2d())
             : before::pixel_list_labels(
                 coords.const_ref(), labeller.size(), labeller.last_frame(), threed);
  };

  for (auto _ : state) {
    af::shared<int> labels = run();
    benchmark::DoNotOptimize(labels.begin());
    benchmark::ClobberMemory();
  }
  set_counters(state, coords.size(), expected);
}

// -------------------------------------------------------------------------
// LabelPixels (hash-table lookup of neighbours)
// -------------------------------------------------------------------------

void BM_LabelPixels(benchmark::State& state, Pattern p, bool after) {
  const int side = state.range(0);
  const int frames = state.range(1);
  LabelPixels label_pixels = make_label_pixels(p, side, frames);
  af::shared<vec3<int> > coords = label_pixels.coords();

  af::shared<int> expected = label_pixels.labels();
  check_same(expected,
             before::label_pixels_labels(coords.const_ref(), label_pixels.size()));
  release_heap();

  auto run = [&]() {
    return after ? label_pixels.labels()
                 : before::label_pixels_labels(coords.const_ref(), label_pixels.size());
  };

  for (auto _ : state) {
    af::shared<int> labels = run();
    benchmark::DoNotOptimize(labels.begin());
    benchmark::ClobberMemory();
  }
  set_counters(state, coords.size(), expected);
}

// Image side and number of frames. A template, since the Benchmark class
// has moved namespace between versions of Google Benchmark.
template <typename B>
B* sizes_3d(B* b) {
  return b->Args({64, 8})->Args({256, 8})->Args({512, 8});
}
template <typename B>
B* sizes_2d(B* b) {
  return b->Args({256, 1})->Args({1024, 1})->Args({2048, 1});
}

/**
 * Frames of 2048x2048 needed for ~16M strong pixels, the size of the
 * synthetic test in https://github.com/dials/dials/issues/3286
 */
int frames_for_16m(Pattern p) {
  switch (p) {
  case Sparse:
    return 192;  // 2% of 192 frames: 16.1M
  case Spots:
    return 72;  // 27 of every 512 pixels, over 72 frames: 15.9M
  case Solid:
    return 4;  // 16.8M
  case Serpentine:
    return 8;  // A little over half of 8 frames: 16.8M
  }
  return 0;
}

int register_benchmarks() {
  const Pattern patterns[] = {Sparse, Spots, Solid, Serpentine};
  const char* when[] = {"Before", "After"};
  for (Pattern p : patterns) {
    for (int after = 0; after < 2; ++after) {
      std::string suffix = std::string(pattern_name(p)) + "/" + when[after];
      sizes_3d(
        benchmark::RegisterBenchmark(
          ("PixelListLabeller3d/" + suffix).c_str(), BM_PixelList, p, true, after))
        ->Args({2048, frames_for_16m(p)})
        ->Unit(benchmark::kMillisecond);
      sizes_2d(
        benchmark::RegisterBenchmark(
          ("PixelListLabeller2d/" + suffix).c_str(), BM_PixelList, p, false, after))
        ->Unit(benchmark::kMillisecond);
      sizes_3d(benchmark::RegisterBenchmark(
                 ("LabelPixels/" + suffix).c_str(), BM_LabelPixels, p, after))
        ->Unit(benchmark::kMillisecond);
    }
  }
  return 0;
}

const int registered = register_benchmarks();

}  // namespace

BENCHMARK_MAIN();
