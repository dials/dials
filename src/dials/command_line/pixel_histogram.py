from __future__ import annotations

import sys

from dials.array_family import flex


def pixel_histogram(spot_filename, maximum=0x10000):
    """Compute a histogram of pixel values up to the given maximum value"""

    spots = flex.reflection_table.from_file(spot_filename)

    if "shoebox" not in spots:
        sys.exit(f"no shoeboxes in {spot_filename}")

    histogram = flex.histogram(
        flex.double(), data_min=0, data_max=maximum, n_slots=maximum
    )

    for j in range(len(spots)):
        box = spots[j]["shoebox"].data
        histogram.update(
            flex.histogram(
                box.as_double().as_1d(), data_min=0, data_max=maximum, n_slots=maximum
            )
        )

    for c, v in zip(histogram.slot_centers(), histogram.slots()):
        print(c, v)


if __name__ == "__main__":
    pixel_histogram(sys.argv[1])
