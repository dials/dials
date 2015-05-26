#!/usr/bin/env python
#
# spot_finder_factory.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def generate_phil_scope():
  from iotbx.phil import parse
  import dials.extensions # import dependency
  from dials.interfaces import SpotFinderThresholdIface

  phil_scope = parse('''

  spotfinder
    .help = "Parameters used in the spot finding algorithm."
  {
    include scope dials.data.lookup.phil_scope
    include scope dials.data.multiprocessing.phil_scope

    write_hot_mask = True
      .type = bool
      .help = "Write the hot mask"

    scan_range = None
      .help = "The range of images to use in finding spots. Number of arguments"
              "must be a factor of two. Specifying \"0 0\" will use all images"
              "by default. The given range follows C conventions"
              "(e.g. j0 <= j < j1)."
              "For sweeps the scan range is interpreted as the literal scan"
              "range. Whereas for imagesets the scan range is interpreted as"
              "the image number in the imageset"
      .type = ints(size=2)
      .multiple = True

    filter
      .help = "Parameters used in the spot finding filter strategy."

    {
      min_spot_size = Auto
        .help = "The minimum number of contiguous pixels for a spot"
                "to be accepted by the filtering algorithm."
        .type = int(value_min=0)

      max_separation = 2
        .help = "The maximum peak-to-centroid separation (in pixels)"
                "for a spot to be accepted by the filtering algorithm."
        .type = float(value_min=0)
        .expert_level = 1

      d_min = None
        .help = "The high resolution limit in Angstrom for a spot to be"
                "accepted by the filtering algorithm."
        .type = float(value_min=0)

      d_max = None
        .help = "The low resolution limit in Angstrom for a spot to be"
                "accepted by the filtering algorithm."
        .type = float(value_min=0)

      background_gradient
        .expert_level=2
      {
        filter = False
          .type = bool
        background_size = 2
          .type = int(value_min=1)
        gradient_cutoff = 4
          .type = float(value_min=0)
      }

      spot_density
        .expert_level=2
      {
        filter = False
          .type = bool
      }

      ice_rings {
        filter = False
          .type = bool
        unit_cell = 4.498,4.498,7.338,90,90,120
          .type = unit_cell
          .help = "The unit cell to generate d_spacings for powder rings."
          .expert_level = 1
        space_group = 194
          .type = space_group
          .help = "The space group used to generate d_spacings for powder rings."
          .expert_level = 1
        width = 0.06
          .type = float(value_min=0.0)
          .help = "The width of an ice ring (in d-spacing)."
          .expert_level = 1
      }

      untrusted_polygon = None
        .multiple = True
        .type = ints(value_min=0)
        .help = "The pixel coordinates (fast, slow) that define the corners "
                "of the untrusted polygon. Spots whose centroids fall within "
                "the bounds of the untrusted polygon will be rejected."

      #untrusted_ellipse = None
      #  .multiple = True
      #  .type = ints(size=4, value_min=0)

      #untrusted_rectangle = None
      #  .multiple = True
      #  .type = ints(size=4, value_min=0)
    }
  }

  ''', process_includes=True)

  main_scope = phil_scope.get_without_substitution("spotfinder")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(SpotFinderThresholdIface.phil_scope())
  return phil_scope

phil_scope = generate_phil_scope()

class FilterRunner(object):
  '''
  A class to run multiple filters in succession.

  '''

  def __init__(self, filters=None):
    '''
    Initialise with a list of filters.

    :param filters: The list of filters

    '''
    if filters is None:
      self.filters=[]
    else:
      self.filters=filters

  def __call__(self, flags, **kwargs):
    '''
    Call the filters one by one.

    :param flags: The input flags
    :returns: The filtered flags

    '''
    flags = self.check_flags(flags, **kwargs)
    for f in self.filters:
      flags = f(flags, **kwargs)
    return flags

  def check_flags(self, flags, predictions=None, observations=None,
                  shoeboxes=None, **kwargs):
    '''
    Check the flags are set, if they're not then create a list
    of Trues equal to the number of items given.

    :param flags: The input flags
    :param predictions: The predictions
    :param observations: The observations
    :param shoeboxes: The shoeboxes
    :return: The filtered flags

    '''
    from scitbx.array_family import flex

    # If flags are not set then create a list of Trues
    if flags == None:
      length = 0
      if predictions:
        length = len(predictions)
      if observations:
        if length > 0:
          assert(length == len(observations))
        else:
          length = len(observations)
      if shoeboxes:
        if length > 0:
          assert(length == len(observations))
        else:
          length = len(shoeboxes)

      # Create an array of flags
      flags = flex.bool(length, True)

    # Return the flags
    return flags


class MinPixelsFilter(object):
  '''
  Filter the reflections by the number of pixels in the shoeboxes.

  '''

  def __init__(self, num, code):
    '''
    Initialise

    :param num: The minimum number of pixels allowed
    :param code: The mask code to use for comparison

    '''
    self.code = code
    self.num = num

  def run(self, flags, observations=None, shoeboxes=None, **kwargs):
    '''
    Run the filtering.

    '''

    # Get the number of mask values matching the code
    count = shoeboxes.count_mask_values(self.code)

    # Return the flags of those > the given number
    return flags.__and__(count >= self.num)

  def __call__(self, flags, **kwargs):
    '''
    Call the filter and print information.

    '''
    from logging import info
    info('Filtering {0} spots by number of pixels'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by number of pixels'.format(
        flags.count(True)))
    return flags


class PeakCentroidDistanceFilter(object):

  def __init__(self, maxd):
    '''
    Initialise

    :param maxd: The maximum distance allowed

    '''
    self.maxd = maxd

  def run(self, flags, observations=None, shoeboxes=None, **kwargs):
    '''
    Run the filtering.

    '''

    # Get the peak locations and the centroids and return the flags of
    # those closer than the min distance
    peak = shoeboxes.peak_coordinates()
    cent = observations.centroids().px_position()
    return flags.__and__((peak - cent).norms() <= self.maxd)

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by peak-centroid distance'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by peak-centroid distance'.format(
        flags.count(True)))
    return flags


class CentroidResolutionFilter(object):

  def __init__(self, d_min, d_max):
    '''
    Initialise

    :param dmin: The maximum resolution
    :param dmax: The minimum resolution

    '''
    if d_min == None:
      self.d_min = 0.0
    else:
      self.d_min = d_min

    if d_max == None:
      self.d_max = 1000.0
    else:
      self.d_max = d_max

  def run(self, flags, sweep=None, observations=None, **kwargs):
    '''
    Run the filtering.

    '''

    # Get all the observation resolutions
    d = observations.resolution(sweep.get_beam(), sweep.get_detector())

    # Return the flags of those in range
    return (flags.__and__(d >= self.d_min)).__and__(d <= self.d_max)

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by resolution'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by resolution'.format(
        flags.count(True)))
    return flags


class PowderRingFilter(object):

  def __init__(self, crystal_symmetry, width=0.06):
    self.crystal_symmetry = crystal_symmetry
    self.width = width

  def run(self, flags, sweep=None, observations=None, **kwargs):
    from cctbx import uctbx

    from dials.array_family import flex
    from dxtbx import imageset
    detector = sweep.get_detector()
    beam = sweep.get_beam()

    ms = self.crystal_symmetry.build_miller_set(
      anomalous_flag=False, d_min=detector.get_max_resolution(beam.get_s0()))
    ms = ms.sort(by_value="resolution")

    miller_indices = flex.miller_index()
    two_thetas_obs = flex.double()
    wavelength = beam.get_wavelength()

    half_width = 0.5 * self.width
    for i, centroid in enumerate(observations):
      if not flags[i]: continue
      x, y = centroid.centroid.px_xy
      d_spacing = detector[centroid.panel].get_resolution_at_pixel(
        beam.get_s0(), (x, y))
      for j, d in enumerate(ms.d_spacings().data()):
        if abs(d - d_spacing) < half_width:
          flags[i] = False
          miller_indices.append(ms.indices()[j])
          two_thetas_obs.append(uctbx.d_star_sq_as_two_theta(
            uctbx.d_as_d_star_sq(d_spacing), wavelength=wavelength, deg=True))

    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by powder rings'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by powder rings'.format(
        flags.count(True)))
    return flags


class polygon(object):
  def __init__(self, vertices):
    assert len(vertices) > 2
    self.vertices = vertices

  def is_inside(self, x, y):
    # http://en.wikipedia.org/wiki/Point_in_polygon
    # http://en.wikipedia.org/wiki/Even-odd_rule
    poly = self.vertices
    num = len(poly)
    i = 0
    j = num - 1
    inside = False
    for i in range(num):
      if  ((poly[i][1] > y) != (poly[j][1] > y)) and \
          (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1]) + poly[i][0]):
        inside = not inside
      j = i
    return inside


class UntrustedPolygonFilter(object):

  def __init__(self, polygons):
    self.polygons = polygons

  def run(self, flags, sweep=None, observations=None, **kwargs):
    for i, centroid in enumerate(observations):
      if not flags[i]: continue
      x, y = centroid.centroid.px_xy
      for poly in self.polygons:
        if poly.is_inside(x, y):
          flags[i] = False
    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by untrusted polygons'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by untrusted polygons'.format(
        flags.count(True)))
    return flags


class BackgroundGradientFilter(object):

  def __init__(self, background_size=2, gradient_cutoff=4):
    self.background_size = background_size
    self.gradient_cutoff = gradient_cutoff

  def run(self, flags, sweep=None, shoeboxes=None, **kwargs):
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    from dials.algorithms.background.simple import Linear2dModeller
    bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
    fg_code = MaskCode.Valid | MaskCode.Foreground
    strong_code = MaskCode.Valid | MaskCode.Strong

    modeller = Linear2dModeller()
    expanded_shoeboxes = flex.shoebox()
    detector = sweep.get_detector()

    zoffset = 0
    if sweep.get_scan() is not None:
      zoffset = sweep.get_scan().get_array_range()[0]

    from libtbx.containers import OrderedDict
    class image_data_cache(object):
      def __init__(self, imageset, size=10):
        self.imageset = imageset
        self.size = size
        self._image_data = OrderedDict()

      def __getitem__(self, i):
        image_data = self._image_data.get(i)
        if image_data is None:
          image_data = self.imageset[i]
          if len(self._image_data) >= self.size:
            # remove the oldest entry in the cache
            del self._image_data[self._image_data.keys()[0]]
          self._image_data[i] = image_data
        return image_data

    cache = image_data_cache(sweep)
    #cache = sweep

    # sort shoeboxes by centroid z
    frame = shoeboxes.centroid_all().position_frame()
    perm = flex.sort_permutation(frame)
    shoeboxes = shoeboxes.select(perm)
    buffer_size = 1
    bg_plus_buffer = self.background_size + buffer_size

    import time
    t0 = time.time()
    for i, shoebox in enumerate(shoeboxes):
      if not flags[perm[i]]:
        continue
      panel = detector[shoebox.panel]
      trusted_range = panel.get_trusted_range()
      max_x, max_y = panel.get_image_size()
      bbox = shoebox.bbox
      x1, x2, y1, y2, z1, z2 = bbox
      # expand the bbox with a background region around the spotfinder shoebox
      # perhaps also should use a buffer zone between the shoebox and the
      # background region
      expanded_bbox = (max(0, x1-bg_plus_buffer),
                       min(max_x, x2+bg_plus_buffer),
                       max(0, y1-bg_plus_buffer),
                       min(max_y, y2+bg_plus_buffer),
                       z1, z2)
      shoebox.bbox = expanded_bbox
    t1 = time.time()
    info("Time expand_shoebox: %s" %(t1-t0))

    rlist = flex.reflection_table()
    rlist['shoebox'] = shoeboxes
    rlist['shoebox'].allocate()
    rlist['panel'] = shoeboxes.panels()
    rlist['bbox'] = shoeboxes.bounding_boxes()

    t0 = time.time()
    rlist.extract_shoeboxes(sweep)
    t1 = time.time()

    shoeboxes = rlist['shoebox']
    shoeboxes.flatten()

    t0 = time.time()
    for i, shoebox in enumerate(shoeboxes):
      if not flags[perm[i]]: continue
      panel = detector[shoebox.panel]
      trusted_range = panel.get_trusted_range()
      max_x, max_y = panel.get_image_size()
      ex1, ex2, ey1, ey2, ez1, ez2 = shoebox.bbox
      data = shoebox.data
      mask = flex.bool(data.accessor(), False)
      for i_y, y in enumerate(range(ey1, ey2)):
        for i_x, x in enumerate(range(ex1, ex2)):
          value = data[0, i_y, i_x]
          if (y >= (ey1+buffer_size) and y < (ey2-buffer_size) and
              x >= (ex1+buffer_size) and x < (ex2-buffer_size)):
            mask[0, i_y, i_x] = False # foreground
          elif (value > trusted_range[0] and value < trusted_range[1]):
            mask[0, i_y, i_x] = True # background

      model = modeller.create(data.as_double(), mask)
      d, a, b = model.params()[:3]
      c = -1

      if abs(a) > self.gradient_cutoff or abs(b) > self.gradient_cutoff:
        flags[perm[i]] = False

      # FIXME should this commented out section be removed?
      #if abs(a) < self.gradient_cutoff and abs(b) < self.gradient_cutoff:
        #flags[i] = False

      #if x2-x1 > 10 or y2-y1 > 10:
        #print a, b, d, flags[perm[i]]
        #bg = flex.double(data.accessor())
        #for x in range(ex2-ex1):
          #for y in range(ey2-ey1):
            #z = a * x + b * y + d
            #bg[0,y,x] = z

        #model = modeller.create(data-bg, mask)
        #d, a, b = model.params()[:3]
        #c = -1

        #bg2 = flex.double(data.accessor())
        #for x in range(ex2-ex1):
          #for y in range(ey2-ey1):
            #z = a * x + b * y + d
            #bg2[0,y,x] = z
        ##print a, b, d

        #from matplotlib import pyplot
        #fig, axes = pyplot.subplots(nrows=1, ncols=5)
        #im0 = axes[0].imshow(data.as_numpy_array()[i_z,:,:], interpolation='none')
        #im1 = axes[1].imshow(mask.as_numpy_array()[i_z,:,:], interpolation='none')
        #im2 = axes[2].imshow(bg.as_numpy_array()[i_z,:,:], interpolation='none')
        #im3 = axes[3].imshow((data-bg).as_numpy_array()[i_z,:,:], interpolation='none')
        #im4 = axes[4].imshow(bg2.as_numpy_array()[i_z,:,:], interpolation='none')
        ##pyplot.colorbar(im0)
        ##pyplot.colorbar(im1)
        ##pyplot.colorbar(im2)
        ##pyplot.colorbar(im3)
        #pyplot.show()

      #from matplotlib import pyplot
      #fig, axes = pyplot.subplots(nrows=1, ncols=2)
      #im0 = axes[0].imshow(data.as_numpy_array()[i_z,:,:], interpolation='none')
      #im1 = axes[1].imshow(bg.as_numpy_array()[i_z,:,:], interpolation='none')
      #pyplot.colorbar(im1)
      #pyplot.show()

    t1 = time.time()
    #print "Time fit_bg: %s" %(t1-t0)

    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by background gradient'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by background gradient'.format(
        flags.count(True)))
    return flags


class SpotDensityFilter(object):

  def __init__(self, nbins=50, gradient_cutoff=0.002):
    self.nbins = nbins
    self.gradient_cutoff = gradient_cutoff

  def run(self, flags, sweep=None, observations=None, **kwargs):
    obs_x, obs_y = observations.centroids().px_position_xy().parts()

    import numpy as np
    H, xedges, yedges = np.histogram2d(
      obs_x.as_numpy_array(), obs_y.as_numpy_array(),bins=self.nbins)

    from scitbx.array_family import flex
    H_flex = flex.double(H.flatten().astype(np.float64))
    n_slots = min(int(flex.max(H_flex)), 30)
    hist = flex.histogram(H_flex, n_slots=n_slots)

    slots = hist.slots()
    cumulative_hist = flex.long(len(slots))
    for i in range(len(slots)):
      cumulative_hist[i] = slots[i]
      if i > 0:
        cumulative_hist[i] += cumulative_hist[i-1]

    cumulative_hist = cumulative_hist.as_double()/flex.max(
      cumulative_hist.as_double())

    cutoff = None
    gradients = flex.double()
    for i in range(len(slots)-1):
      x1 = cumulative_hist[i]
      x2 = cumulative_hist[i+1]
      g = (x2 - x1)/hist.slot_width()
      gradients.append(g)
      if (cutoff is None and  i > 0 and
          g < self.gradient_cutoff and gradients[i-1] < self.gradient_cutoff):
        cutoff = hist.slot_centers()[i-1]-0.5*hist.slot_width()

    H_flex = flex.double(np.ascontiguousarray(H))
    isel = (H_flex > cutoff).iselection()
    sel = np.column_stack(np.where(H > cutoff))
    for (ix, iy) in sel:
      flags.set_selected(
        ((obs_x > xedges[ix]) & (obs_x < xedges[ix+1]) &
         (obs_y > yedges[iy]) & (obs_y < yedges[iy+1])), False)

    if 0:
      from matplotlib import pyplot
      fig, ax1 = pyplot.subplots()
      extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
      plot1 = ax1.imshow(H, extent=extent, interpolation="nearest")
      pyplot.xlim((0, pyplot.xlim()[1]))
      pyplot.ylim((0, pyplot.ylim()[1]))
      pyplot.gca().invert_yaxis()
      cbar1 = pyplot.colorbar(plot1)
      pyplot.axes().set_aspect('equal')
      pyplot.show()

      fig, ax1 = pyplot.subplots()
      ax2 = ax1.twinx()
      ax1.scatter(hist.slot_centers()-0.5*hist.slot_width(), cumulative_hist)
      ax1.set_ylim(0, 1)
      ax2.plot(hist.slot_centers()[:-1]-0.5*hist.slot_width(), gradients)
      ymin, ymax = pyplot.ylim()
      pyplot.vlines(cutoff, ymin, ymax, color='r')
      pyplot.show()

      H2 = H.copy()
      if cutoff is not None:
        H2[np.where(H2 >= cutoff)] = 0
      fig, ax1 = pyplot.subplots()
      plot1 = ax1.pcolormesh(xedges, yedges, H2)
      pyplot.xlim((0, pyplot.xlim()[1]))
      pyplot.ylim((0, pyplot.ylim()[1]))
      pyplot.gca().invert_yaxis()
      cbar1 = pyplot.colorbar(plot1)
      pyplot.axes().set_aspect('equal')
      pyplot.show()

    return flags

  def __call__(self, flags, **kwargs):
    ''' Call the filter and print information. '''
    from logging import info
    info('Filtering {0} spots by spot density'.format(
        flags.count(True)))
    flags = self.run(flags, **kwargs)
    info('Filtered {0} spots by spot density'.format(
        flags.count(True)))
    return flags


class SpotFinderFactory(object):
  '''
  Factory class to create spot finders

  '''

  @staticmethod
  def from_parameters(params=None):
    '''
    Given a set of parameters, construct the spot finder

    :param params: The input parameters
    :returns: The spot finder instance

    '''
    from dials.algorithms.peak_finding.spot_finder import SpotFinder
    from libtbx.phil import parse

    if params is None:
      params = phil_scope.fetch(source=parse("")).extract()

    # Read in the lookup files
    mask = SpotFinderFactory.load_image(params.spotfinder.lookup.mask)
    params.spotfinder.lookup.mask = mask

    # Configure the algorithm and wrap it up
    find_spots = SpotFinderFactory.configure_algorithm(params)
    filter_spots = SpotFinderFactory.configure_filter(params)
    return SpotFinder(
      find_spots=find_spots,
      filter_spots=filter_spots,
      scan_range=params.spotfinder.scan_range,
      write_hot_mask=params.spotfinder.write_hot_mask)

  @staticmethod
  def configure_algorithm(params):
    '''
    Given a set of parameters, construct the spot finder

    :param params: The input parameters
    :returns: The spot finder instance

    '''
    from dials.algorithms.peak_finding.spot_finder import ExtractSpots

    # Create the threshold strategy
    threshold = SpotFinderFactory.configure_threshold(params)

    # Setup the spot finder
    return ExtractSpots(threshold_image=threshold,
                        mask=params.spotfinder.lookup.mask,
                        mp_method=params.spotfinder.mp.method,
                        nproc=params.spotfinder.mp.nproc)

  @staticmethod
  def configure_threshold(params):
    '''
    Get the threshold strategy

    :param params: The input parameters
    :return: The threshold algorithm

    '''
    from dials.interfaces import SpotFinderThresholdIface
    Algorithm = SpotFinderThresholdIface.extension(
      params.spotfinder.threshold.algorithm)
    return Algorithm(params)

  @staticmethod
  def configure_filter(params):
    '''
    Get the filter strategy.

    :param params: The input parameters
    :return: The filter algorithm

    '''
    from dials.algorithms import shoebox
    from cctbx import crystal

    # Initialise an empty list of filters
    filters = []

    # Add a min number of pixels filter
    if params.spotfinder.filter.min_spot_size is not None:
      filters.append(MinPixelsFilter(
          params.spotfinder.filter.min_spot_size,
          shoebox.MaskCode.Valid))

    # Add a peak-centroid distance filter
    if params.spotfinder.filter.max_separation is not None:
      filters.append(PeakCentroidDistanceFilter(
        params.spotfinder.filter.max_separation))

    # Add a centroid resolution filter
    if (params.spotfinder.filter.d_min is not None or
        params.spotfinder.filter.d_max is not None):
      filters.append(CentroidResolutionFilter(
          params.spotfinder.filter.d_min,
          params.spotfinder.filter.d_max))

    if params.spotfinder.filter.ice_rings.filter:
      crystal_symmetry = crystal.symmetry(
        unit_cell=params.spotfinder.filter.ice_rings.unit_cell,
        space_group=params.spotfinder.filter.ice_rings.space_group.group())
      filters.append(
        PowderRingFilter(crystal_symmetry,
                         width=params.spotfinder.filter.ice_rings.width))

    if len(params.spotfinder.filter.untrusted_polygon):
      polygons = []
      for vertices in params.spotfinder.filter.untrusted_polygon:
        if vertices is not None:
          assert len(vertices) % 2 == 0
          vertices = [vertices[i*2:i*2+2] for i in range(len(vertices)//2)]
          polygons.append(polygon(vertices))
      if len(polygons):
        filters.append(UntrustedPolygonFilter(polygons))

    if params.spotfinder.filter.background_gradient.filter:
      bg_filter_params = params.spotfinder.filter.background_gradient
      filters.append(BackgroundGradientFilter(
        background_size=bg_filter_params.background_size,
        gradient_cutoff=bg_filter_params.gradient_cutoff))

    if params.spotfinder.filter.spot_density.filter:
      filters.append(SpotDensityFilter())

    # Return the filter runner with the list of filters
    return FilterRunner(filters)

  @staticmethod
  def load_image(filename_or_data):
    '''
    Given a filename, load an image. If the data is already loaded, return it.

    :param filename_or_data: The input filename (or data)
    :return: The image or None

    '''
    import cPickle as pickle

    # If no filename is set then return None
    if not filename_or_data:
      return None

    # If it's already loaded, return early
    if type(filename_or_data) is tuple:
      return filename_or_data

    # Read the image and return the image data
    image = pickle.load(open(filename_or_data))
    if not isinstance(image, tuple):
      image = (image,)
    return image
