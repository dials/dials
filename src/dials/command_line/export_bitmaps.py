from __future__ import annotations

import pathlib
import sys

import PIL.Image

import iotbx.phil

from dials.util import export_bitmaps, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """

Export raw diffraction image files as bitmap images, optionally exporting
images from intermediate spot-finding steps (local mean and variance maps,
or sigma_b, sigma_s or threshold-filtered images). Appearance of the images
can be altered via the brightness and colour_scheme parameters, and optionally
binning of pixels can be used to reduce image sizes.

Examples::

  dials.export_bitmaps image.cbf

  dials.export_bitmaps models.expt

  dials.export_bitmaps image.cbf display=variance colour_scheme=inverse_greyscale
"""

phil_scope = iotbx.phil.parse(
    """
binning = 1
  .type = int(value_min=1)
brightness = 100
  .type = float(value_min=0.0)
colour_scheme = *greyscale rainbow heatmap inverse_greyscale
  .type = choice
projection = lab *image
  .type = choice
padding = 4
  .type = int(value_min=0)
imageset_index = None
  .type = int
  .multiple = True
  .help = "The index/indices from an imageset to export. The first image of "
          "the set is 1."
  .expert_level=2
display = *image mean variance dispersion sigma_b \
          sigma_s threshold global_threshold
  .type = choice
nsigma_b = 6
  .type = float(value_min=0)
nsigma_s = 3
  .type = float(value_min=0)
global_threshold = 0
  .type = float(value_min=0)
kernel_size = 3,3
  .type = ints(size=2, value_min=1)
min_local = 2
  .type = int
gain = 1
  .type = float(value_min=0)
saturation = 0
  .type = int
show_mask = False
  .type = bool
resolution_rings {
  show = False
    .type = bool
  number = 5
    .type = int(value_min=1)
  d_spacings = None
    .type = floats(value_min=0)
  fontsize = 30
    .type = int
    .optional = True
  fill = red
    .type = str
    .help = "Colour of the resolution rings and labels"
}
ice_rings {
  show = False
    .type = bool
  fontsize = None
    .type = int
    .optional = True
  unit_cell = 4.498,4.498,7.338,90,90,120
    .type = unit_cell
    .help = "The unit cell to generate d_spacings for powder rings."
    .expert_level = 1
  space_group = 194
    .type = space_group
    .help = "The space group used to generate d_spacings for powder rings."
    .expert_level = 1
  fill = blue
    .type = str
    .help = "Colour of the ice rings and labels"
}
png {
  compress_level = 1
    .type = int(value_min=0, value_max=9)
    .help = "ZLIB compression level, a number between 0 and 9: 1 gives best "
            "speed, 9 gives best compression, 0 gives no compression at all."
}
jpeg {
  quality = 75
    .type = int(value_min=1, value_max=95)
    .help = "The image quality, on a scale from 1 (worst) to 95 (best)"
}

include scope dials.util.options.format_phil_scope

output {
  prefix = "image"
    .type = str
  directory = None
    .type = path
  file = None
    .type = str
    .help = "Full name of the output file. Overrides 'prefix' and the default "
            "file extension. Only makes sense if a single file is written."
  format = jpeg *png tiff
    .type = choice
}""",
    process_includes=True,
)

colour_schemes = {"greyscale": 0, "rainbow": 1, "heatmap": 2, "inverse_greyscale": 3}


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.export_bitmaps [options] models.expt | image.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        check_format=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    imagesets = experiments.imagesets()

    for imageset in imagesets:
        imageset_as_bitmaps(imageset, params)


def imageset_as_bitmaps(imageset, params):
    if params.output.directory is None:
        params.output.directory = "."
    output_dir = pathlib.Path(params.output.directory)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    output_files = []

    if (
        scan := imageset.get_scan()
    ) is not None and not scan.is_still():  # and not images:
        start, end = scan.get_image_range()
    else:
        start, end = 1, len(imageset)
    images = list(range(start, end + 1))
    if params.imageset_index:
        selected_images = []
        for idx in params.imageset_index:
            if idx > len(images):
                sys.exit(
                    f"Bad value for imageset_index: {idx}. The imageset has length {len(images)}; allowable values for imageset_index are the range 1 to {len(images)} inclusive."
                )
            selected_images.append(images[idx - 1])
        images = selected_images

    if params.output.file and len(images) != 1:
        sys.exit("output.file can only be specified if a single image is exported")

    for i_img, flex_img in enumerate(
        export_bitmaps.imageset_as_flex_image(
            imageset,
            images,
            brightness=params.brightness,
            binning=params.binning,
            projection=export_bitmaps.Projection(params.projection),
            saturation=params.saturation,
            show_mask=params.show_mask,
            display=export_bitmaps.Display(params.display),
            colour_scheme=export_bitmaps.ColourScheme[params.colour_scheme.upper()],
            gain=params.gain,
            nsigma_b=params.nsigma_b,
            global_threshold=params.global_threshold,
            min_local=params.min_local,
            kernel_size=params.kernel_size,
        )
    ):
        pil_img = PIL.Image.frombytes(
            "RGB", (flex_img.ex_size2(), flex_img.ex_size1()), flex_img.as_bytes()
        )
        if params.resolution_rings.show:
            export_bitmaps.draw_resolution_rings(
                imageset,
                pil_img,
                flex_img,
                n_rings=params.resolution_rings.number,
                spacings=params.resolution_rings.d_spacings,
                fill=params.resolution_rings.fill,
                fontsize=params.resolution_rings.fontsize,
                binning=params.binning,
            )
        if params.ice_rings.show:
            export_bitmaps.draw_ice_rings(
                imageset,
                pil_img,
                flex_img,
                unit_cell=params.ice_rings.unit_cell,
                space_group=params.ice_rings.space_group.group(),
                fill=params.ice_rings.fill,
                fontsize=params.ice_rings.fontsize,
                binning=params.binning,
            )
        if params.output.file:
            path = output_dir / params.output.file
        else:
            path = (
                output_dir
                / f"{params.output.prefix}{images[i_img]:0{params.padding}}.{params.output.format}"
            )

        print(f"Exporting {path}")
        output_files.append(path)

        pil_img.save(
            path,
            format=params.output.format,
            compress_level=params.png.compress_level,
            quality=params.jpeg.quality,
        )

    return output_files


if __name__ == "__main__":
    run()
