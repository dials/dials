/*
 * compressed_slice_store.h
 *
 *  Copyright (C) 2026 Diamond Light Source
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_COMPRESSED_SLICE_STORE_H
#define DIALS_ALGORITHMS_INTEGRATION_COMPRESSED_SLICE_STORE_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>
#include <blosc2.h>
#include <lz4.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * The bytes of the data, mask and background of a decoded shoebox
   * @param npixels The number of pixels in the shoebox
   */
  inline std::size_t decoded_shoebox_bytes(std::size_t npixels) {
    return npixels * (2 * sizeof(float) + sizeof(uint8_t));
  }

  /**
   * Choose how many bytes of compressed shoeboxes to decode and process at
   * once on one frame. Usually all the shoeboxes that finish on a frame are
   * processed together. On a few frames many more finish than usual (on the
   * last frame of a scan, every shoebox that runs off the end), and decoding
   * all of those at once would set the peak memory. So the limit is set above
   * a typical frame (1.5 times the median over frames where any shoeboxes
   * finish), and only those outliers are processed in batches.
   * @param finishing The decoded bytes of the shoeboxes finishing on each frame
   * @returns The most bytes to decode at once
   */
  inline std::size_t finishing_batch_bytes(std::vector<std::size_t> finishing) {
    const std::size_t min_bytes = std::size_t(64) << 20;
    finishing.erase(std::remove(finishing.begin(), finishing.end(), 0),
                    finishing.end());
    if (finishing.empty()) {
      return min_bytes;
    }
    std::vector<std::size_t>::iterator median =
      finishing.begin() + finishing.size() / 2;
    std::nth_element(finishing.begin(), median, finishing.end());
    return std::max(min_bytes, *median + *median / 2);
  }

  /**
   * Compressed storage for the z-slices of shoeboxes that are still being read.
   *
   * Each reflection has a stream of slice records, appended one at a time and
   * read back all at once. A record holds the float data, bitshuffled and then
   * compressed with LZ4, and the mask, run-length encoded. Bitshuffle puts each
   * bit position of the pixels together, so the bits that are always zero for
   * small counts become long runs that LZ4 removes. Nothing depends on the data
   * being integers.
   *
   * The streams are kept in fixed-size chunks from a pool owned by the store.
   * Freed chunks are reused rather than returned to malloc, so the many small,
   * short-lived allocations don't fragment the heap.
   */
  class CompressedSliceStore {
  public:
    /** Bytes of each pooled chunk, including its link to the next chunk */
    static const std::size_t chunk_size = 256;

    /** Chunks per slab allocated from the heap */
    static const std::size_t slab_chunks = 4096;

    /** Fixed bytes per reflection in the job */
    static const std::size_t bytes_per_reflection = 12;

    /**
     * Assumed bytes per pixel of a stored slice (data and mask), for memory
     * estimates. Integer counts from photon-counting detectors take about
     * 0.5 B/px and quantised CCD data about 1.3 B/px. Data that don't compress
     * are stored raw, at up to 4 B/px plus the mask.
     */
    static double assumed_bytes_per_pixel() {
      return 2.0;
    }

    /** @param nrefl The number of reflections that may be stored */
    explicit CompressedSliceStore(std::size_t nrefl)
        : streams_(nrefl), free_(npos), nchunks_(0), nused_(0), max_used_(0) {
      blosc2_init();
    }

    /**
     * Append a slice to the stream of a reflection
     * @param index The reflection index
     * @param data The slice data
     * @param mask The slice mask
     * @param n The number of pixels in the slice
     */
    void append(std::size_t index,
                const float* data,
                const uint8_t* mask,
                std::size_t n) {
      DIALS_ASSERT(index < streams_.size());
      DIALS_ASSERT(n > 0);
      std::size_t nbytes = n * sizeof(float);
      DIALS_ASSERT(nbytes < (std::size_t)std::numeric_limits<int>::max() / 2);
      reserve(n);

      // Bitshuffle, then compress. Keep the raw floats if that doesn't help.
      int32_t status = blosc2_bitshuffle(sizeof(float), nbytes, data, shuffled_.data());
      DIALS_ASSERT(status == (int32_t)nbytes);
      uint8_t* out = record_.data() + 4;
      int csize = LZ4_compress_default(reinterpret_cast<const char*>(shuffled_.data()),
                                       reinterpret_cast<char*>(out),
                                       (int)nbytes,
                                       (int)nbytes - 1);
      uint32_t header;
      if (csize > 0) {
        header = (uint32_t)csize << 1;
      } else {
        std::memcpy(out, data, nbytes);
        csize = nbytes;
        header = ((uint32_t)nbytes << 1) | 1;
      }
      std::memcpy(record_.data(), &header, 4);
      std::size_t size = 4 + csize;

      // Run-length encode the mask as (value, varint length) pairs
      for (std::size_t i = 0; i < n;) {
        std::size_t j = run_end(mask, i, n);
        record_[size++] = mask[i];
        for (std::size_t run = j - i; true; run >>= 7) {
          if (run < 0x80) {
            record_[size++] = (uint8_t)run;
            break;
          }
          record_[size++] = (uint8_t)(run | 0x80);
        }
        i = j;
      }
      write(streams_[index], record_.data(), size);
    }

    /**
     * Decode all the slices stored for a reflection, and free its storage.
     * @param index The reflection index
     * @param data The output data, of nslices * n values
     * @param mask The output mask, of nslices * n values
     * @param n The number of pixels in each slice
     * @param nslices The number of slices to decode
     */
    void extract(std::size_t index,
                 float* data,
                 uint8_t* mask,
                 std::size_t n,
                 std::size_t nslices) {
      DIALS_ASSERT(index < streams_.size());
      Stream& stream = streams_[index];
      if (nslices == 0) {
        DIALS_ASSERT(stream.head == npos);
        return;
      }
      std::size_t nbytes = n * sizeof(float);
      reserve(n);

      // Gather the stream into one buffer, freeing the chunks as we go
      std::size_t size = 0;
      for (uint32_t c = stream.head; c != npos;) {
        Chunk& chunk = get(c);
        std::size_t used = (c == stream.tail) ? stream.tail_used : payload_size;
        if (gather_.size() < size + used) {
          gather_.resize(2 * (size + used));
        }
        std::memcpy(gather_.data() + size, chunk.payload, used);
        size += used;
        uint32_t next = chunk.next;
        release(c);
        c = next;
      }
      stream = Stream();

      // Decode each slice record in turn
      const uint8_t* p = gather_.data();
      const uint8_t* end = p + size;
      for (std::size_t z = 0; z < nslices; ++z) {
        float* sdata = data + z * n;
        uint8_t* smask = mask + z * n;
        DIALS_ASSERT(p + 4 <= end);
        uint32_t header;
        std::memcpy(&header, p, 4);
        p += 4;
        std::size_t csize = header >> 1;
        DIALS_ASSERT(p + csize <= end);
        if (header & 1) {
          DIALS_ASSERT(csize == nbytes);
          std::memcpy(sdata, p, nbytes);
        } else {
          int dsize = LZ4_decompress_safe(reinterpret_cast<const char*>(p),
                                          reinterpret_cast<char*>(shuffled_.data()),
                                          (int)csize,
                                          (int)nbytes);
          DIALS_ASSERT(dsize == (int)nbytes);
          int32_t status =
            blosc2_bitunshuffle(sizeof(float), nbytes, shuffled_.data(), sdata);
          DIALS_ASSERT(status == (int32_t)nbytes);
        }
        p += csize;
        for (std::size_t i = 0; i < n;) {
          DIALS_ASSERT(p < end);
          uint8_t value = *p++;
          std::size_t run = 0;
          for (int shift = 0; true; shift += 7) {
            DIALS_ASSERT(p < end);
            uint8_t b = *p++;
            run |= (std::size_t)(b & 0x7f) << shift;
            if (!(b & 0x80)) {
              break;
            }
          }
          DIALS_ASSERT(run > 0 && i + run <= n);
          std::memset(smask + i, value, run);
          i += run;
        }
      }
      DIALS_ASSERT(p == end);
    }

    /** @returns The bytes in chunks currently holding data */
    std::size_t bytes_used() const {
      return nused_ * chunk_size;
    }

    /** @returns The most bytes ever in use */
    std::size_t max_bytes_used() const {
      return max_used_ * chunk_size;
    }

    /** @returns The bytes allocated from the heap for chunks */
    std::size_t bytes_allocated() const {
      return nchunks_ * chunk_size;
    }

  private:
    static const uint32_t npos = 0xffffffff;
    static const std::size_t payload_size = chunk_size - sizeof(uint32_t);

    struct Chunk {
      uint32_t next;
      uint8_t payload[payload_size];
    };

    struct Stream {
      uint32_t head;
      uint32_t tail;
      uint32_t tail_used;
      Stream() : head(npos), tail(npos), tail_used(0) {}
    };

    /** @returns The end of the run of equal bytes starting at i */
    static std::size_t run_end(const uint8_t* mask, std::size_t i, std::size_t n) {
      // Compare 8 bytes at a time while we can
      const uint64_t pattern = 0x0101010101010101ULL * mask[i];
      std::size_t j = i + 1;
      for (; j + 8 <= n; j += 8) {
        uint64_t word;
        std::memcpy(&word, mask + j, 8);
        uint64_t diff = word ^ pattern;
        if (diff != 0) {
          // The first differing byte, in memory order
          const uint8_t* bytes = reinterpret_cast<const uint8_t*>(&diff);
          while (*bytes == 0) {
            ++bytes;
          }
          return j + (bytes - reinterpret_cast<const uint8_t*>(&diff));
        }
      }
      while (j < n && mask[j] == mask[i]) {
        ++j;
      }
      return j;
    }

    Chunk& get(uint32_t c) {
      return slabs_[c / slab_chunks][c % slab_chunks];
    }

    uint32_t acquire() {
      if (free_ == npos) {
        DIALS_ASSERT(nchunks_ + slab_chunks < npos);
        slabs_.push_back(std::unique_ptr<Chunk[]>(new Chunk[slab_chunks]));
        for (std::size_t i = slab_chunks; i > 0; --i) {
          uint32_t c = nchunks_ + i - 1;
          get(c).next = free_;
          free_ = c;
        }
        nchunks_ += slab_chunks;
      }
      uint32_t c = free_;
      Chunk& chunk = get(c);
      free_ = chunk.next;
      chunk.next = npos;
      if (++nused_ > max_used_) {
        max_used_ = nused_;
      }
      return c;
    }

    void release(uint32_t c) {
      get(c).next = free_;
      free_ = c;
      DIALS_ASSERT(nused_ > 0);
      --nused_;
    }

    void write(Stream& stream, const uint8_t* src, std::size_t size) {
      if (stream.head == npos) {
        stream.head = stream.tail = acquire();
        stream.tail_used = 0;
      }
      while (size > 0) {
        if (stream.tail_used == payload_size) {
          uint32_t c = acquire();
          get(stream.tail).next = c;
          stream.tail = c;
          stream.tail_used = 0;
        }
        std::size_t n = std::min(size, payload_size - stream.tail_used);
        std::memcpy(get(stream.tail).payload + stream.tail_used, src, n);
        stream.tail_used += n;
        src += n;
        size -= n;
      }
    }

    /** Make sure the scratch buffers can hold a slice of n pixels */
    void reserve(std::size_t n) {
      std::size_t nbytes = n * sizeof(float);
      // Header, raw data, and up to two RLE bytes per pixel plus slack
      std::size_t record_bytes = 4 + nbytes + 2 * n + 16;
      if (shuffled_.size() < nbytes) {
        shuffled_.resize(nbytes);
      }
      if (record_.size() < record_bytes) {
        record_.resize(record_bytes);
      }
    }

    std::vector<Stream> streams_;
    std::vector<std::unique_ptr<Chunk[]> > slabs_;
    uint32_t free_;
    std::size_t nchunks_;
    std::size_t nused_;
    std::size_t max_used_;
    std::vector<uint8_t> shuffled_;
    std::vector<uint8_t> record_;
    std::vector<uint8_t> gather_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_COMPRESSED_SLICE_STORE_H
