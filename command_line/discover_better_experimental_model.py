from __future__ import division

import math
from libtbx.phil import command_line
import iotbx.phil
from scitbx import matrix
from cctbx.array_family import flex
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks, flatten_reflections
from dials.algorithms.indexing.indexer \
     import indexer_base, filter_reflections_by_scan_range

master_phil_scope = iotbx.phil.parse("""
nproc = Auto
  .type = int(value_min=1)
plot_search_scope = False
  .type = bool
scan_range = None
  .help = "The range of images to use in indexing. Number of arguments"
    "must be a factor of two. Specifying \"0 0\" will use all images"
    "by default. The given range follows C conventions"
    "(e.g. j0 <= j < j1)."
  .type = ints(size=2)
  .multiple = True
mm_search_scope = 4.0
  .help = "Global radius of origin offset search."
  .type = float(value_min=0)
wide_search_binning = 5
  .help = "Modify the coarseness of the wide grid search for the beam centre."
  .type = float(value_min=0)
output = optimized_datablock.json
  .type = path
""")

master_params = master_phil_scope.fetch().extract()


class better_experimental_model_discovery(object):
  def __init__(self, imagesets, spot_lists, solution_lists,
               amax_lists, horizon_phil, wide_search_binning=1):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

  def optimize_origin_offset_local_scope(self):
    """Local scope: find the optimal origin-offset closest to the current overall detector position
       (local minimum, simple minimization)"""

    from libtbx.test_utils import approx_equal
    from rstbx.indexing_api import dps_extended

    detector = self.imagesets[0].get_detector()
    beam = self.imagesets[0].get_beam()
    goniometer = self.imagesets[0].get_goniometer()
    s0 = matrix.col(beam.get_s0())
    # construct two vectors that are perpendicular to the beam.  Gives a basis for refining beam
    if goniometer is None:
      axis = matrix.col((1,0,0))
    else:
      axis = matrix.col(goniometer.get_rotation_axis())
    beamr0 = s0.cross(axis).normalize()
    beamr1 = beamr0.cross(s0).normalize()
    beamr2 = beamr1.cross(s0).normalize()
    #beamr1 = matrix.col(detector[0].get_fast_axis())
    #beamr2 = matrix.col(detector[0].get_slow_axis())
    #beamr0 = beamr1.cross(beamr2).normalize()

    assert approx_equal(s0.dot(beamr1), 0.)
    assert approx_equal(s0.dot(beamr2), 0.)
    assert approx_equal(beamr2.dot(beamr1), 0.)
    # so the orthonormal vectors are self.S0_vector, beamr1 and beamr2

    if self.horizon_phil.indexing.mm_search_scope:
      scope = self.horizon_phil.indexing.mm_search_scope
      plot_px_sz = self.imagesets[0].get_detector()[0].get_pixel_size()[0]
      plot_px_sz *= self.wide_search_binning
      grid = max(1,int(scope/plot_px_sz))
      widegrid = 2 * grid + 1
      scores = flex.double()
      for y in xrange(-grid,grid+1):
        for x in xrange(-grid,grid+1):
          new_origin_offset = x*plot_px_sz*beamr1 + y*plot_px_sz*beamr2
          score = 0
          for i in range(len(self.imagesets)):
            score += self.get_origin_offset_score(new_origin_offset,
                                                  self.solution_lists[i],
                                                  self.amax_lists[i],
                                                  self.spot_lists[i],
                                                  self.imagesets[i])
          scores.append(score)

      plot_max = flex.max(scores)
      idx_max = flex.max_index(scores)

      def igrid(x): return x - (widegrid//2)
      idxs = [igrid(i)*plot_px_sz for i in xrange(widegrid)]

      wide_search_offset = (idxs[idx_max%widegrid])*beamr1 + (idxs[idx_max//widegrid])*beamr2

    else:
      wide_search_offset = None

    # DO A SIMPLEX MINIMIZATION
    from scitbx.simplex import simplex_opt
    class test_simplex_method(object):
      def __init__(selfOO, wide_search_offset=None):
        selfOO.starting_simplex=[]
        selfOO.n = 2
        selfOO.wide_search_offset = wide_search_offset
        for ii in range(selfOO.n+1):
          selfOO.starting_simplex.append(flex.random_double(selfOO.n))
        selfOO.optimizer = simplex_opt( dimension=selfOO.n,
                                      matrix  = selfOO.starting_simplex,
                                      evaluator = selfOO,
                                      tolerance=1e-7)
        selfOO.x = selfOO.optimizer.get_solution()
        selfOO.offset = selfOO.x[0]*0.2*beamr1 + selfOO.x[1]*0.2*beamr2
        if selfOO.wide_search_offset is not None:
          selfOO.offset += selfOO.wide_search_offset

      def target(selfOO, vector):
        trial_origin_offset = vector[0]*0.2*beamr1 + vector[1]*0.2*beamr2
        if selfOO.wide_search_offset is not None:
          trial_origin_offset += selfOO.wide_search_offset
        target = 0
        for i in range(len(self.imagesets)):
          target -= self.get_origin_offset_score(trial_origin_offset,
                                                 self.solution_lists[i],
                                                 self.amax_lists[i],
                                                 self.spot_lists[i],
                                                 self.imagesets[i])
        return target

    MIN = test_simplex_method(wide_search_offset=wide_search_offset)
    new_offset = MIN.offset

    if self.horizon_phil.indexing.plot_search_scope:
      scope = self.horizon_phil.indexing.mm_search_scope
      plot_px_sz = self.imagesets[0].get_detector()[0].get_pixel_size()[0]
      grid = max(1,int(scope/plot_px_sz))
      scores = flex.double()
      for y in xrange(-grid,grid+1):
        for x in xrange(-grid,grid+1):
          new_origin_offset = x*plot_px_sz*beamr1 + y*plot_px_sz*beamr2
          score = 0
          for i in range(len(self.imagesets)):
            score += self.get_origin_offset_score(new_origin_offset,
                                                  self.solution_lists[i],
                                                  self.amax_lists[i],
                                                  self.spot_lists[i],
                                                  self.imagesets[i])
          scores.append(score)

      def show_plot(widegrid,excursi):
        excursi.reshape(flex.grid(widegrid, widegrid))
        plot_max = flex.max(excursi)
        idx_max = flex.max_index(excursi)

        def igrid(x): return x - (widegrid//2)
        idxs = [igrid(i)*plot_px_sz for i in xrange(widegrid)]

        from matplotlib import pyplot as plt
        plt.figure()
        CS = plt.contour([igrid(i)*plot_px_sz for i in xrange(widegrid)],
                         [igrid(i)*plot_px_sz for i in xrange(widegrid)], excursi.as_numpy_array())
        plt.clabel(CS, inline=1, fontsize=10, fmt="%6.3f")
        plt.title("Wide scope search for detector origin offset")
        plt.scatter([0.0],[0.0],color='g',marker='o')
        plt.scatter([new_offset[0]] , [new_offset[1]],color='r',marker='*')
        plt.scatter([idxs[idx_max%widegrid]] , [idxs[idx_max//widegrid]],color='k',marker='s')
        plt.axes().set_aspect("equal")
        plt.xlabel("offset (mm) along beamr1 vector")
        plt.ylabel("offset (mm) along beamr2 vector")
        plt.show()

        #changing value
        trial_origin_offset =  (idxs[idx_max%widegrid])*beamr1 + (idxs[idx_max//widegrid])*beamr2
        return trial_origin_offset

      show_plot(widegrid = 2 * grid + 1, excursi = scores)

    return dps_extended.get_new_detector(self.imagesets[0].get_detector(), new_offset)

  def get_origin_offset_score(self, trial_origin_offset, solutions, amax, spots_mm, imageset):
    from rstbx.indexing_api import lattice # import dependency
    from rstbx.indexing_api import dps_extended
    trial_detector = dps_extended.get_new_detector(imageset.get_detector(), trial_origin_offset)

    from dials.algorithms.indexing.indexer import indexer_base
    indexer_base.map_centroids_to_reciprocal_space(
      spots_mm, trial_detector, imageset.get_beam(), imageset.get_goniometer())

    return self.sum_score_detail(spots_mm['rlp'], solutions, amax=amax)

  def sum_score_detail(self, reciprocal_space_vectors, solutions, granularity=None, amax=None):
    """Evaluates the probability that the trial value of ( S0_vector | origin_offset ) is correct,
       given the current estimate and the observations.  The trial value comes through the
       reciprocal space vectors, and the current estimate comes through the short list of
       DPS solutions. Actual return value is a sum of NH terms, one for each DPS solution, each ranging
       from -1.0 to 1.0"""
    import cmath
    from rstbx.dps_core import Direction, Directional_FFT
    nh = min ( solutions.size(), 20) # extended API
    #solutions = self.getSolutions() #extended API
    sum_score = 0.0
    for t in xrange(nh):
      #if t!=unique:continue
      dfft = Directional_FFT(
        angle=Direction(solutions[t]), xyzdata=reciprocal_space_vectors,
        granularity=5.0, amax=amax, # extended API XXX These values have to come from somewhere!
        F0_cutoff = 11)
      kval = dfft.kval();
      kmax = dfft.kmax();
      #kval_cutoff = self.raw_spot_input.size()/4.0; # deprecate record
      kval_cutoff = reciprocal_space_vectors.size()/4.0; # deprecate record
      if ( kval > kval_cutoff ):
        ff=dfft.fft_result;
        kbeam = ((-dfft.pmin)/dfft.delta_p) + 0.5;
        Tkmax = cmath.phase(ff[kmax]);
        backmax = math.cos(Tkmax+(2*math.pi*kmax*kbeam/(2*ff.size()-1)) );
        ### Here it should be possible to calculate a gradient.
        ### Then minimize with respect to two coordinates.  Use lbfgs?  Have second derivatives?
        ### can I do something local to model the cosine wave?
        ### direction of wave travel.  Period. phase.
        sum_score += backmax;
      #if t == unique:
      #  print t, kmax, dfft.pmin, dfft.delta_p, Tkmax,(2*math.pi*kmax*kbeam/(2*ff.size()-1))
    return sum_score


def run_dps(args):
  imageset, spots, params = args

  detector = imageset.get_detector()
  beam = imageset.get_beam()
  goniometer = imageset.get_goniometer()
  scan = imageset.get_scan()

  spots_mm = indexer_base.map_spots_pixel_to_mm_rad(
    spots=spots, detector=detector, scan=scan)

  #from dials.algorithms.indexing.indexer import Indexer
  #spots_mm = Indexer._map_spots_pixel_to_mm_rad(
    #spots=spots, detector=detector, scan=scan)

  # derive a max_cell from mm spots
  # derive a grid sampling from spots

  from rstbx.indexing_api.lattice import DPS_primitive_lattice
  # max_cell: max possible cell in Angstroms; set to None, determine from data
  # recommended_grid_sampling_rad: grid sampling in radians; guess for now

  DPS = DPS_primitive_lattice(max_cell = None,
                              recommended_grid_sampling_rad = None,
                              horizon_phil = params)
  from scitbx import matrix
  DPS.S0_vector = matrix.col(beam.get_s0())
  DPS.inv_wave = 1./beam.get_wavelength()
  if goniometer is None:
    DPS.axis = matrix.col((1,0,0))
  else:
    DPS.axis = matrix.col(goniometer.get_rotation_axis())
  DPS.set_detector(detector)

  # transform input into what Nick needs
  # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

  data = flex.vec3_double()
  for spot in spots_mm:
    data.append((spot['xyzobs.mm.value'][0],
                 spot['xyzobs.mm.value'][1],
                 spot['xyzobs.mm.value'][2]*180./math.pi))

  #from matplotlib import pyplot as plt
  #plt.plot([spot.centroid_position[0] for spot in spots_mm] , [spot.centroid_position[1] for spot in spots_mm], 'ro')
  #plt.show()

  print "Running DPS"

  DPS.index(raw_spot_input=data,
            panel_addresses=flex.int([s['panel'] for s in spots_mm]))
  print "Found %i solutions" %len(DPS.getSolutions()),
  print "with max unit cell %7.2f Angstroms."%(DPS.amax)
  return dict(solutions=flex.vec3_double([s.dvec for s in DPS.getSolutions()]),amax=DPS.amax)


def discover_better_experimental_model(imagesets, spot_lists, params,
                                       nproc=1, wide_search_binning=1):
  assert len(imagesets) == len(spot_lists)
  assert len(imagesets) > 0
  # XXX should check that all the detector and beam objects are the same
  from dials.algorithms.indexing.indexer import indexer_base
  spot_lists_mm = [
    indexer_base.map_spots_pixel_to_mm_rad(
      spots, imageset.get_detector(), imageset.get_scan())
    for spots, imageset in zip(spot_lists, imagesets)]

  args = [(imageset, spots, params)
          for imageset, spots in zip(imagesets, spot_lists_mm)]

  from libtbx import easy_mp
  results = easy_mp.parallel_map(
    func=run_dps,
    iterable=args,
    processes=nproc,
    method="multiprocessing",
    preserve_order=True,
    asynchronous=True,
    preserve_exception_message=True)
  solution_lists = [r["solutions"] for r in results]
  amax_list = [r["amax"] for r in results]
  assert len(solution_lists) > 0

  detector = imagesets[0].get_detector()
  beam = imagesets[0].get_beam()

  # perform calculation
  if params.indexing.improve_local_scope == "origin_offset":
    discoverer = better_experimental_model_discovery(
      imagesets, spot_lists_mm, solution_lists, amax_list, params,
      wide_search_binning=wide_search_binning)
    new_detector = discoverer.optimize_origin_offset_local_scope()
    old_beam_centre = detector.get_ray_intersection(beam.get_s0())[1]
    new_beam_centre = new_detector.get_ray_intersection(beam.get_s0())[1]
    print "Old beam centre: %.2f mm, %.2f mm" %old_beam_centre
    print "New beam centre: %.2f mm, %.2f mm" %new_beam_centre
    print "Shift: %.2f mm, %.2f mm" %(
      matrix.col(old_beam_centre)-matrix.col(new_beam_centre)).elems
    return new_detector, beam
  elif params.indexing.improve_local_scope=="S0_vector":
    raise NotImplementedError()
    #new_S0_vector = DPS.optimize_S0_local_scope()
    #import copy
    #new_beam = copy.copy(beam)
    #new_beam.set_s0(new_S0_vector)
    #return detector, new_beam


def run(args):
  if len(args) == 0:
    from libtbx.utils import Usage
    import libtbx.load_env
    from cStringIO import StringIO
    usage_message = """\
%s datablock.json strong.pickle [options]

Parameters:
""" %libtbx.env.dispatcher_name
    s = StringIO()
    master_phil_scope.show(out=s)
    usage_message += s.getvalue()
    raise Usage(usage_message)

  parser = OptionParser(
    phil=master_phil_scope,
    read_datablocks=True,
    read_reflections=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)
  if len(datablocks) == 0:
    print "No DataBlock could be constructed"
    return
  imagesets = []
  for datablock in datablocks:
    imagesets.extend(datablock.extract_imagesets())

  assert len(imagesets) > 0
  assert len(reflections) == len(imagesets)

  if params.scan_range is not None and len(params.scan_range) > 0:
    reflections = [
      filter_reflections_by_scan_range(refl, params.scan_range)
      for refl in reflections]

  from rstbx.phil.phil_preferences import indexing_api_defs
  import iotbx.phil
  hardcoded_phil = iotbx.phil.parse(
    input_string=indexing_api_defs).extract()
  # for development, we want an exhaustive plot of beam probability map:
  hardcoded_phil.indexing.plot_search_scope = params.plot_search_scope
  hardcoded_phil.indexing.mm_search_scope = params.mm_search_scope

  new_detector, new_beam = discover_better_experimental_model(
    imagesets, reflections, hardcoded_phil, nproc=params.nproc,
    wide_search_binning=params.wide_search_binning)
  for imageset in imagesets:
    imageset.set_detector(new_detector)
    imageset.set_beam(new_beam)
  from dxtbx.serialize import dump
  dump.datablock(datablock, params.output)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
