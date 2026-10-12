# Benchmarks

C++ microbenchmarks using [Google Benchmark](https://github.com/google/benchmark).
They are not built by default. To build them, configure with:

    cmake -DDIALS_BUILD_BENCHMARKS=ON ...
    ninja bench_connected_components
    ./dials/benchmarks/bench_connected_components

An installed Google Benchmark (e.g. `conda install benchmark`) is used if
found; otherwise it is downloaded at configure time.

Build in Release mode (`-DCMAKE_BUILD_TYPE=Release`) for meaningful numbers.

## bench_connected_components

Connected component labelling as used by `dials.find_spots`
(`PixelListLabeller::labels_3d`/`labels_2d`) and `LabelPixels`. Each benchmark
has a `Before` variant, which runs a copy of the old
`boost::adjacency_list` + `boost::connected_components` code, and an `After`
variant, which runs the current union-find (`DisjointSets`). Before and
after are checked to give identical labels before timing.

The inputs are contrived stacks of `side × side` images over a number of frames
(the two benchmark arguments):

| Pattern      | Shape                                                  |
| ------------ | ------------------------------------------------------ |
| `Sparse`     | ~2% of pixels at random; mostly isolated pixels        |
| `Spots`      | 3×3×3 blobs on an 8-pixel grid; many small spots       |
| `Solid`      | every pixel set; one spot, with as many edges as possible |
| `Serpentine` | one pixel-wide snake through each frame; one long, thin spot |

`PixelListLabeller3d` also has a ~16M strong pixel case for each pattern
(2048×2048 frames, with the number of frames chosen per pattern), matching the
size of the synthetic test in
[#3286](https://github.com/dials/dials/issues/3286). Select these with
`--benchmark_filter='PixelListLabeller3d/.*/2048/'`. The `Before` variants
need several GB of RAM at this size.

Use e.g. `--benchmark_filter=PixelListLabeller3d` to select a subset.
