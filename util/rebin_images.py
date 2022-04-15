from __future__ import annotations

import binascii
import os
import random

from dxtbx.ext import compress


def gz_open(filename, mode):
    import gzip

    return gzip.GzipFile(filename, mode)


def split_counts(image, split):
    from scitbx.array_family import flex

    new_images = [flex.int(flex.grid(image.focus()), 0) for k in range(split)]

    negative = image.as_1d() < 0
    positive = image.as_1d() > 0

    for new_image in new_images:
        new_image.as_1d().set_selected(negative, image.as_1d())

    for p in positive.iselection():
        counts = image[p]
        for j in range(counts):
            new_images[random.randint(0, split - 1)][p] += 1

    return new_images


def merge_counts(images):
    from scitbx.array_family import flex

    image = flex.int(flex.grid(images[0].focus()), 0)
    negative = images[0].as_1d() < 0
    for i in images:
        image += i
    image.as_1d().set_selected(negative, images[0].as_1d())
    return image


def read_image(in_image):
    from dxtbx import load

    assert os.path.exists(in_image)

    start_tag = binascii.unhexlify("0c1a04d5")

    data = gz_open(in_image, "rb").read()
    data_offset = data.find(start_tag)
    cbf_header = data[:data_offset]
    pixel_values = load(in_image).get_raw_data()

    return pixel_values, cbf_header


def write_image(out_image, pixel_values, header, nn=1):
    assert not os.path.exists(out_image)
    start_tag = binascii.unhexlify("0c1a04d5")

    compressed = compress(pixel_values)

    fixed_header = ""
    header = header.decode()
    for record in header.split("\n")[:-1]:
        if "X-Binary-Size:" in record:
            fixed_header += f"X-Binary-Size: {len(compressed)}\r\n"
        elif "Content-MD5" in record:
            pass
        elif "Count_cutoff" in record:
            cutoff = int(record.split()[2]) * nn
            fixed_header += "# Count_cutoff %d counts\n" % cutoff
        else:
            fixed_header += f"{record}\n"

    tailer = "\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n"

    gz_open(out_image, "wb").write(
        fixed_header.encode() + start_tag + compressed + tailer.encode()
    )


def main(in_images, out_images):
    assert len(in_images) == len(out_images)
    n = len(in_images)

    for i in in_images:
        assert os.path.exists(i)
    for o in out_images:
        assert not os.path.exists(o)

    in_image_data = []
    in_image_headers = []

    for i in in_images:
        print(f"Reading {i}")
        pixel, header = read_image(i)
        in_image_data.append(pixel)
        in_image_headers.append(header)

    sum_image = merge_counts(in_image_data)
    rebin_images = split_counts(sum_image, n)

    for o, pixel, header in zip(out_images, rebin_images, in_image_headers):
        print(f"Writing {o}")
        write_image(o, pixel, header)


def main_sum(in_images, out_image):
    for i in in_images:
        assert os.path.exists(i)
    assert not os.path.exists(out_image)

    in_image_data = []
    in_image_headers = []

    for i in in_images:
        print(f"Reading {i}")
        pixel, header = read_image(i)
        in_image_data.append(pixel)
        in_image_headers.append(header)

    sum_image = merge_counts(in_image_data)

    print(f"Writing {out_image}")
    write_image(out_image, sum_image, in_image_headers[0], nn=len(in_images))
