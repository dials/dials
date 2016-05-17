from __future__ import division

import iotbx.phil


from dials.util.options import flatten_datablocks

help_message = '''

Export raw diffraction image files as bitmap images, optionally exporting
images from intermediate spot-finding steps (local mean and variance maps,
or sigma_b, sigma_s or threshold-filtered images). Appearance of the images
can be altered via the brightness and colour_scheme parameters, and optionally
binning of pixels can be used to reduce image sizes.

Examples::

  dials.export_bitmaps image.cbf

  dials.export_bitmaps datablock.json

  dials.export_bitmaps image.cbf display=variance colour_scheme=inverse_greyscale

'''

phil_scope = iotbx.phil.parse("""
binning = 1
  .type = int(value_min=1)
brightness = 100
  .type = float(value_min=0.0)
colour_scheme = *greyscale rainbow heatmap inverse_greyscale
  .type = choice
format = jpeg *png tiff
  .type = choice
prefix = "image"
  .type = str
output_dir = None
  .type = path
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
""")

colour_schemes = {
  'greyscale': 0,
  'rainbow': 1,
  'heatmap': 2,
  'inverse_greyscale': 3
}

def run(args):
  import os
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util.options import OptionParser
  usage = "%s [options] datablock.json | image.cbf" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  datablocks = flatten_datablocks(params.input.datablock)
  if len(datablocks) == 0:
    parser.print_help()
    exit(0)

  imagesets = datablocks[0].extract_imagesets()

  brightness = params.brightness / 100
  vendortype = "made up"

  # check that binning is a power of 2
  binning = params.binning
  if not (binning > 0 and ((binning & (binning - 1)) == 0)):
    raise Sorry("binning must be a power of 2")

  output_dir = params.output_dir
  if output_dir is None:
    output_dir = "."
  elif not os.path.exists(output_dir):
    os.makedirs(output_dir)

  from rstbx.slip_viewer.tile_generation \
       import _get_flex_image, _get_flex_image_multipanel

  for imageset in imagesets:
    detector = imageset.get_detector()
    panel = detector[0]
    # XXX is this inclusive or exclusive?
    saturation = panel.get_trusted_range()[1]
    if params.saturation:
      saturation = params.saturation
    for i_image, image in enumerate(imageset):

      if len(detector) == 1:
        image = [image]

      trange = [p.get_trusted_range() for p in detector]
      mask = []
      mask = imageset.get_mask(i_image)
      if mask is None:
        mask = [p.get_trusted_range_mask(im) for im, p in zip(image, detector)]

      image = image_filter(image, mask, display=params.display, gain_value=params.gain,
                           nsigma_b=params.nsigma_b,
                           nsigma_s=params.nsigma_s,
                           global_threshold=params.global_threshold,
                           min_local=params.min_local,
                           kernel_size=params.kernel_size)

      if len(detector) > 1:
        # FIXME This doesn't work properly, as flex_image.size2() is incorrect
        # also binning doesn't work
        assert binning == 1
        flex_image = _get_flex_image_multipanel(
          brightness=brightness,
          panels=detector,
          raw_data=image)
      else:
        flex_image = _get_flex_image(
          brightness=brightness,
          data=image[0],
          binning=binning,
          saturation=saturation,
          vendortype=vendortype)

      flex_image.setWindow(0, 0, 1)
      flex_image.adjust(color_scheme=colour_schemes.get(params.colour_scheme))

      # now export as a bitmap
      flex_image.prep_string()
      import Image
      # XXX is size//binning safe here?
      pil_img = Image.fromstring(
        'RGB', (flex_image.size2()//binning,
                flex_image.size1()//binning),
        flex_image.export_string)

      path = os.path.join(
        output_dir, params.prefix + ("%04d" % i_image) + '.' + params.format)

      print "Exporting %s" %path
      with open(path, 'wb') as tmp_stream:
        pil_img.save(tmp_stream, format=params.format)

def image_filter(raw_data, mask, display,
                 gain_value, nsigma_b, nsigma_s, global_threshold,
                 min_local, kernel_size):

  from dials.algorithms.image.threshold import KabschDebug
  from dials.array_family import flex

  if display == 'image':
    return raw_data

  assert gain_value > 0
  gain_map = [flex.double(raw_data[i].accessor(), gain_value)
              for i in range(len(raw_data))]

  kabsch_debug_list = []
  for i_panel in range(len(raw_data)):
    kabsch_debug_list.append(
      KabschDebug(
        raw_data[i_panel].as_double(), mask[i_panel], gain_map[i_panel],
        kernel_size, nsigma_b, nsigma_s, global_threshold, min_local))

  if display == 'mean':
    display_data = [kabsch.mean() for kabsch in kabsch_debug_list]
  elif display == 'variance':
    display_data = [kabsch.variance() for kabsch in kabsch_debug_list]
  elif display == 'dispersion':
    display_data = [
      kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
  elif display == 'sigma_b':
    cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
    display_data = [kabsch.cv_mask() for kabsch in kabsch_debug_list]
    display_data = [mask.as_1d().as_double() for mask in display_data]
    for i, mask in enumerate(display_data):
      mask.reshape(cv[i].accessor())
  elif display == 'sigma_s':
    cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
    display_data = [kabsch.value_mask() for kabsch in kabsch_debug_list]
    display_data = [mask.as_1d().as_double() for mask in display_data]
    for i, mask in enumerate(display_data):
      mask.reshape(cv[i].accessor())
  elif display == 'global_threshold':
    cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
    display_data = [kabsch.global_mask() for kabsch in kabsch_debug_list]
    display_data = [mask.as_1d().as_double() for mask in display_data]
    for i, mask in enumerate(display_data):
      mask.reshape(cv[i].accessor())
  elif display == 'threshold':
    cv = [kabsch.coefficient_of_variation() for kabsch in kabsch_debug_list]
    display_data = [kabsch.final_mask() for kabsch in kabsch_debug_list]
    display_data = [mask.as_1d().as_double() for mask in display_data]
    for i, mask in enumerate(display_data):
      mask.reshape(cv[i].accessor())

  return display_data


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
