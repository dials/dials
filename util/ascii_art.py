from __future__ import absolute_import, division, print_function

import six


def spot_counts_per_image_plot(reflections, **kwargs):
    if len(reflections) == 0:
        return "\n"

    x, y, z = reflections["xyzobs.px.value"].parts()
    return flex_histogram(z, **kwargs)


def histogram_from_json(filename):
    from dials.array_family import flex
    import json

    with open(filename, "r") as fh:
        z = json.load(fh)
    return flex_histogram(flex.double(z))


def flex_histogram(z, char="*", width=60, height=10):
    from dials.array_family import flex

    assert isinstance(char, six.string_types)
    assert len(char) == 1

    # import json
    # with open('list.json', 'w') as fh:
    #   json.dump(sorted(list(z)), fh)

    min_z = flex.min(z)
    max_z = flex.max(z)

    epsilon = flex.double(n % 2 for n in range(len(z)))
    epsilon = (epsilon / 2) - 0.25
    z += epsilon

    # image numbers to display on x-axis label
    xmin = int(round(min_z + 1e-16))
    xmax = int(round(max_z))

    # estimate the total number of images
    image_count = xmax - xmin + 1
    if image_count <= 1:
        return "%i spots found on 1 image" % len(z)

    # determine histogram width
    width = min(image_count, width)
    z_step = (max_z - min_z) / width

    # bin all spots
    counts = flex.double()
    seen_spots = 0
    for i in range(1, width):
        z_bound = min_z + (i * z_step)
        spots = (z < z_bound).count(True)
        counts.append(spots - seen_spots)
        seen_spots = spots
    spots = (z >= z_bound).count(True)
    counts.append(spots)

    max_count = flex.max(counts)
    total_counts = flex.sum(counts)
    assert total_counts == len(
        z
    ), "Only found %d out of %d reflections for histogram" % (total_counts, len(z))
    counts *= height / max_count
    counts = counts.iround()

    rows = []
    rows.append(
        "%i spots found on %i images (max %i / bin)"
        % (total_counts, image_count, max_count)
    )

    for i in range(height, 0, -1):
        row = []
        for j, c in enumerate(counts):
            if c > (i - 1):
                row.append(char)
            else:
                row.append(" ")
        rows.append("".join(row))

    padding = width - len(str(xmin)) - len(str(xmax))
    rows.append(
        "%i%s%i"
        % (
            xmin,
            (" " if padding < 7 else "image").center(padding) if padding > 0 else "",
            xmax,
        )
    )
    return "\n".join(rows)


if __name__ == "__main__":
    import sys
    from libtbx import easy_pickle

    for arg in sys.argv[1:]:
        print(arg)
        print(spot_counts_per_image_plot(easy_pickle.load(arg)))
        print()
