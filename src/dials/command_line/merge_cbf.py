from __future__ import annotations

import binascii
import sys

import iotbx.phil
from dxtbx.ext import compress, uncompress

import dials.util

help_message = """

This program can be used to merge a given number of consecutive cbf files into
a smaller number of images For example, running dials.merge_cbf on a experiments
with 100 images, using the default value of merge_n_images=2, will output 50
summed images, with every consecutive pair of images being summed into a single
output image. Currently only cbf format images are supported as input.

Examples::

  dials.merge_cbf image_*.cbf

  dials.merge_cbf image_*.cbf merge_n_images=10
"""

phil_scope = iotbx.phil.parse(
    """\
merge_n_images = 2
  .type = int(value_min=1)
  .help = "Number of input images to average into a single output image"

get_raw_data_from_imageset = True
  .type = bool
  .help = "By default the raw data is read via the imageset. This limits use"
          "to single panel detectors where the format class does not make"
          "modifications to the array size in the file. Set this option to"
          "false in order to bypass the imageset and read the data as-is from"
          "the CBF file"
  .expert_level = 2

output {
  image_prefix = sum_
    .type = path
}
""",
    process_includes=True,
)


def get_raw_data_from_file(imageset, i):
    """Read cbf directly to access the raw data array rather than
    through the imageset, in order to work for multi-panel detectors and other
    situations where the format class modifies the raw array"""
    file_name = imageset.get_image_identifier(i)
    with open(file_name, "rb") as cbf:
        data = cbf.read()
    start_tag = binascii.unhexlify("0c1a04d5")
    data_offset = data.find(start_tag) + 4
    cbf_header = data[: data_offset - 4].decode("latin-1")
    fast = slow = length = 0
    for record in cbf_header.split("\n"):
        if "X-Binary-Size-Fastest-Dimension" in record:
            fast = int(record.split()[-1])
        elif "X-Binary-Size-Second-Dimension" in record:
            slow = int(record.split()[-1])
        elif "X-Binary-Size:" in record:
            length = int(record.split()[-1])
    values = uncompress(
        packed=data[data_offset : data_offset + length], fast=fast, slow=slow
    )
    return (values,)


def merge_cbf(imageset, n_images, out_prefix="sum_", get_raw_data_from_imageset=True):
    from dxtbx.format.FormatCBF import FormatCBF

    assert issubclass(
        imageset.get_format_class(), FormatCBF
    ), "Only CBF format images supported"

    assert len(imageset) >= n_images

    n_output_images = len(imageset) // n_images

    n_digits = len(str(n_output_images))

    for i_out in range(n_output_images):
        data_out = None

        for j in range(n_images):

            i_in = (i_out * n_images) + j

            if get_raw_data_from_imageset:
                data_in = imageset.get_raw_data(i_in)
            else:
                data_in = get_raw_data_from_file(imageset, i_in)

            assert len(data_in) == 1
            data_in = data_in[0]
            if data_out is None:
                data_out = data_in
            else:
                # FIXME only add pixels to this which are > 0; image pixels < 0
                # are meaningful and should be preserved;
                # Achieved by setting -ve values here to 0 before +=
                # This assumes that -ve values are constant over all images
                data_special = data_in < 0
                data_in.set_selected(data_special, 0)
                data_out += data_in

        out_image = "{prefix}{number:0{digits}d}.cbf".format(
            prefix=out_prefix, number=i_out + 1, digits=n_digits
        )

        start_tag = binascii.unhexlify("0c1a04d5")

        with open(imageset.get_path(i_out * n_images), "rb") as fh:
            data = fh.read()
        data_offset = data.find(start_tag)
        cbf_header = data[:data_offset].decode("latin-1")

        new_header = []
        compressed = compress(data_out)

        old_size = 0

        for record in cbf_header.split("\n")[:-1]:
            rsplit = record.split(" ")
            if "X-Binary-Size:" in record:
                old_size = int(record.split()[-1])
                new_header.append(f"X-Binary-Size: {len(compressed)}\r\n")
            elif "Content-MD5" in record:
                pass
            elif len(rsplit) > 3 and rsplit[1] in {
                "Exposure_time",
                "Angle_increment",
                "Exposure_period",
                "Count_cutoff",
                "Phi_increment",
                "Omega_increment",
                "Chi_increment",
            }:

                if rsplit[1] == "Count_cutoff":  # needs to be an integer
                    new_header.append(
                        "%s\n"
                        % " ".join(
                            rsplit[:2]
                            + ["%d" % (n_images * int(rsplit[2]))]
                            + rsplit[3:]
                        )
                    )
                else:
                    new_header.append(
                        "%s\n"
                        % " ".join(
                            rsplit[:2]
                            + ["%f" % (n_images * float(rsplit[2]))]
                            + rsplit[3:]
                        )
                    )

            else:
                new_header.append(f"{record}\n")

        loop_lines = [
            n for n, record in enumerate(new_header) if record.startswith("loop_")
        ]
        multiply_fields = {
            "_diffrn_scan_axis.angle_range",
            "_diffrn_scan_axis.angle_increment",
            "_diffrn_scan_axis.displacement_range",
            "_diffrn_scan_axis.displacement_increment",
            "_diffrn_scan_frame.integration_time",
            "_diffrn_scan_frame.exposure_time",
            "_array_intensities.overload",
        }
        for loop_start in loop_lines:
            n = loop_start
            modifiers = []
            while True:
                n = n + 1
                line = new_header[n].strip()
                if line in {"", ";"}:  # end of loop
                    break
                elif line.startswith("_"):  # loop header
                    if line in multiply_fields:
                        modifiers.append(n_images)
                    else:
                        modifiers.append(None)
                elif any(modifiers):  # loop body
                    # NOTE: This can break when fields are modified in loops with
                    #   'Strings with spaces, as they are seen as multiple columns, or with'
                    #   _multiple _columns _defined _on _same _line _they _are _seen _as _one _column
                    new_line = [
                        element
                        if modifier is None
                        else f"{float(element) * modifier:f}"
                        for modifier, element in zip(modifiers, line.split())
                    ]
                    new_header[n] = f"{' '.join(new_line)}\r\n"

        tailer = data[data_offset + 4 + old_size :]

        with open(out_image, "wb") as f:
            f.write("".join(new_header).encode("latin-1"))
            f.write(start_tag)
            f.write(compressed)
            f.write(tailer)
        print(f"{out_image} written")


@dials.util.show_mail_handle_errors()
def run(args=None):
    from dials.util.options import ArgumentParser, flatten_experiments

    usage = "dials.merge_cbf [options] image_*.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args, show_diff_phil=True, return_unhandled=True
    )

    n_images = params.merge_n_images
    out_prefix = params.output.image_prefix
    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        return

    if len(experiments) > 1:
        sys.exit("Only one experiment can be processed at a time")
    else:
        imagesets = experiments.imagesets()
        assert len(imagesets) == 1
        imageset = imagesets[0]

    merge_cbf(
        imageset,
        n_images,
        out_prefix=out_prefix,
        get_raw_data_from_imageset=params.get_raw_data_from_imageset,
    )


if __name__ == "__main__":
    run()
