#!/usr/bin/env python

from __future__ import division

from dxtbx.model.goniometer import goniometer_factory
import cPickle as pickle
from dials.array_family import flex
from scitbx import matrix
import sys
import tools

# Create the phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    profile_model = 'profile_model.phil'
      .type = str
      .help = "The profile parameters output filename"

    reflections = 'integrated.pickle'
      .type = str
      .help = "The integrated output filename"
  }

  sampling
    .expert_level = 1
  {

    reflections_per_degree = 50
      .help = "The number of predicted reflections per degree of the sweep "
              "to integrate."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input "
              "prediction list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of predictions to integrate."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    integrate_all_reflections = True
      .help = "Override reflections_per_degree and integrate all predicted"
              "reflections."
      .type = bool

  }

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

''', process_includes=True)


class Script(object):
  """ The integration program. """

  def __init__(self):
    """Initialise the script."""
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the parser
    self.parser = OptionParser(
      usage="you are on your own.",
      phil=phil_scope,
      epilog="you are on your own.",
      read_experiments=True)

  def list_possible_reflections(self, spacegroup, unit_cell, dmin, dmax):
    # List the possible miller indices in the minimum resolution range
    # and subtract those outside the maximum resolution range
    from dials.algorithms.spot_prediction import IndexGenerator

    # TODO: Does this make sense? Pass spacegroup.type() and then filter (again?)?
    possible_indices_p1 = set(IndexGenerator(unit_cell, spacegroup.type(), dmin).to_array())
    lowres_indices_p1 =   set(IndexGenerator(unit_cell, spacegroup.type(), dmax + 1e-8).to_array())
    # [dmin, dmax] are inclusive ranges, thus reduce inner sphere radius by a small value (1e-8)

    possible_indices_p1 = possible_indices_p1 - lowres_indices_p1
    print "%5d unique reflections possible in P 1" % len(possible_indices_p1)

    # filter systematic absent miller indices
    systematic_present_indices = [n for n in possible_indices_p1 if not spacegroup.is_sys_absent(n)]
    print "%5d of these reflections are not systematically absent in %s" %\
          (len(systematic_present_indices), spacegroup.type().lookup_symbol())

    return systematic_present_indices


  def add_seen_multiplicity(self, hkl_seen_multiplicity, hkl_observations):
    # Count the miller indices in hkl_observations and add to the counts recorded in hkl_seen_multiplicity.
    overall_multiplicity = hkl_seen_multiplicity.copy()
    for hkl in hkl_observations:
      if overall_multiplicity.has_key(hkl):
        overall_multiplicity[hkl] += 1
      else:
        overall_multiplicity[hkl] = 1
    return overall_multiplicity


  def StrategyEvaluator(self, experiments, evaluation_function_factory, dmin, dmax):
    # TODO: Don't take experiments, rather take one set of experiment components and a strategy equivalent list

    def calculate_miller_rings(detector):
      import time
      t = time.clock()
      mrings = {}
      for scan_axis in ['phi', 'omega']:
        mrings[scan_axis] = {}
        for other_axis in range(0, 360, 10):
          gonio = goniometer_factory.make_kappa_goniometer(alpha=54.7356, # fixed magic angle
                                                           kappa=0,       # fixed 0 kappa
                                                           phi=0 if scan_axis == "phi" else other_axis,
                                                           omega=0 if scan_axis == "omega" else other_axis,
                                                           direction="-y",
                                                           scan_axis=scan_axis)
          ring = tools.determine_miller_ring_sectors(detector, gonio, s0, possible_hkl_flex, crystal_A)
          mrings[scan_axis][other_axis] = ring

      print time.clock() - t
      print list(mrings.iterkeys())
      print [ list(x.iterkeys()) for x in mrings.itervalues() ]

      t = time.clock()
      sc = scorings
      all_scorings = { '%s:%s:%s' % (scan, str(x1), str(x2)): tools.geometric_scoring([hkl_to_id[hkl] for hkl in mrings[scan][x1][x2]], sc)['total']
                       for scan in mrings.iterkeys() for x1 in mrings[scan].iterkeys() for x2 in range(0,36) }
#     print list(all_scorings.itervalues())
      print len(list(all_scorings.itervalues())), min(all_scorings.itervalues()), max(all_scorings.itervalues()), sum(all_scorings.itervalues()) / len(list(all_scorings.itervalues()))
      print time.clock() - t

      sys.exit(0)

    def calculate_observations(detector, goniometer, oscillation):
      calculate_miller_rings(detector)
      import time
      t = time.clock()
      r = tools.determine_miller_ring_sectors(detector, goniometer, s0, possible_hkl_flex, crystal_A)

      sc = scorings
      print map(lambda(x): tools.geometric_scoring([hkl_to_id[hkl] for hkl in x], sc)['total'], r)
      sc = tools.geometric_scoring([hkl_to_id[hkl] for hkl in r[4]], scorings)['scorings']
      sc = tools.geometric_scoring([hkl_to_id[hkl] for hkl in r[5]], sc)['scorings']
      sc = tools.geometric_scoring([hkl_to_id[hkl] for hkl in r[6]], sc)['scorings']
      print map(lambda(x): tools.geometric_scoring([hkl_to_id[hkl] for hkl in x], sc)['total'], r)
      print time.clock() -t


      print map(len, r)
      sys.exit(0)
      # TODO: To test 2theta rotation code use detector.py detector factory
      # TODO: and create two different two_theta detectors and try to rotate one onto the other
      # ----: Reflection only counted once within oscillation, even when multiple times in diffraction condition
      #   ----: Is it?! Unsure.
      #     SOLVED: half-true. Not problematic here.

    def export_mtz(observed_hkls, experiment, filename):
      if experiment.goniometer:
        axis = experiment.goniometer.get_rotation_axis()
      else:
        axis = 0.0, 0.0, 0.0
      s0 = experiment.beam.get_s0()
      wavelength = experiment.beam.get_wavelength()

      from scitbx import matrix

      panel = experiment.detector[0]
      pixel_size = panel.get_pixel_size()
      cb_op_to_ref = experiment.crystal.get_space_group().info(
        ).change_of_basis_op_to_reference_setting()

      experiment.crystal = experiment.crystal.change_basis(cb_op_to_ref)

      from iotbx import mtz
      from scitbx.array_family import flex
      import itertools

      m = mtz.object()
      m.set_title('from dials.scratch.mg.strategy_i19')
      m.set_space_group_info(experiment.crystal.get_space_group().info())

      nrefcount = sum(observed_hkls.itervalues())
      nref = max(observed_hkls.itervalues())

      for batch in range(1, nref+1):
        o = m.add_batch().set_num(batch).set_nbsetid(1).set_ncryst(1)
        o.set_time1(0.0).set_time2(0.0).set_title('Batch %d' % batch)
        o.set_ndet(1).set_theta(flex.float((0.0, 0.0))).set_lbmflg(0)
        o.set_alambd(wavelength).set_delamb(0.0).set_delcor(0.0)
        o.set_divhd(0.0).set_divvd(0.0)
        o.set_so(flex.float(s0)).set_source(flex.float((0, 0, -1)))
        o.set_bbfac(0.0).set_bscale(1.0)
        o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)
        _unit_cell = experiment.crystal.get_unit_cell()
        _U = experiment.crystal.get_U()

        o.set_cell(flex.float(_unit_cell.parameters()))
        o.set_lbcell(flex.int((-1, -1, -1, -1, -1, -1)))
        o.set_umat(flex.float(_U.transpose().elems))
        mosaic = experiment.crystal.get_mosaicity()
        o.set_crydat(flex.float([mosaic, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        o.set_lcrflg(0)
        o.set_datum(flex.float((0.0, 0.0, 0.0)))

        # detector size, distance
        o.set_detlm(flex.float([0.0, panel.get_image_size()[0],
                                0.0, panel.get_image_size()[1],
                                0, 0, 0, 0]))
        o.set_dx(flex.float([panel.get_distance(), 0.0]))

        # goniometer axes and names, and scan axis number, and number of axes, missets
        o.set_e1(flex.float(axis))
        o.set_e2(flex.float((0.0, 0.0, 0.0)))
        o.set_e3(flex.float((0.0, 0.0, 0.0)))
        o.set_gonlab(flex.std_string(('AXIS', '', '')))
        o.set_jsaxs(1)
        o.set_ngonax(1)
        o.set_phixyz(flex.float((0.0, 0.0, 0.0, 0.0, 0.0, 0.0)))

        phi_start, phi_range = 0.0, 0.0
        o.set_phistt(phi_start)
        o.set_phirange(phi_range)
        o.set_phiend(phi_start + phi_range)
        o.set_scanax(flex.float(axis))

        # number of misorientation angles
        o.set_misflg(0)

        # crystal axis closest to rotation axis (why do I want this?)
        o.set_jumpax(0)

        # type of data - 1; 2D, 2; 3D, 3; Laue
        o.set_ldtype(2)

      # now create the actual data structures - first keep a track of the columns
      # H K L M/ISYM BATCH I SIGI IPR SIGIPR FRACTIONCALC XDET YDET ROT WIDTH
      # LP MPART FLAG BGPKRATIOS

      from cctbx.array_family import flex as cflex # implicit import

      # now go for it and make an MTZ file...
      x = m.add_crystal('XTAL', 'DIALS', unit_cell.parameters())
      d = x.add_dataset('FROMDIALS', wavelength)

      # now add column information...
      type_table = {'IPR': 'J', 'BGPKRATIOS': 'R', 'WIDTH': 'R', 'I': 'J',
                    'H': 'H', 'K': 'H', 'MPART': 'I', 'L': 'H', 'BATCH': 'B',
                    'M_ISYM': 'Y', 'SIGI': 'Q', 'FLAG': 'I', 'XDET': 'R', 'LP': 'R',
                    'YDET': 'R', 'SIGIPR': 'Q', 'FRACTIONCALC': 'R', 'ROT': 'R'}

      m.adjust_column_array_sizes(nrefcount)
      m.set_n_reflections(nrefcount)

      # assign H, K, L, M_ISYM space
      for column in 'H', 'K', 'L', 'M_ISYM':
        d.add_column(column, type_table[column]).set_values(flex.float(nrefcount, 0.0))

      batchnums = ( _ for (x, n) in observed_hkls.iteritems() for _ in range(1, n+1) )
      d.add_column('BATCH', type_table['BATCH']).set_values(flex.float(batchnums))
      d.add_column('FRACTIONCALC', type_table['FRACTIONCALC']).set_values(flex.float(nrefcount, 3.0))

      m.replace_original_index_miller_indices(cb_op_to_ref.apply(
        cflex.miller_index([ _ for (x, n) in observed_hkls.iteritems() for _ in itertools.repeat(x, n) ])
        ))

      m.write(filename)

      return m

    def evaluate_concrete_strategy(hkls, evaluation_function):
      print "%5d reflections fall on detector during sweep" % len(hkls)
#      print "%5d unique reflections fall on detector during sweep (without symmetry relations)" % len(set(list(hkls)))
#      for r in detectable_rays.rows():
#        hkl = r['miller_index']
#        if (abs(hkl[0]) == 11) and (abs(hkl[1]) == 5) and (abs(hkl[2]) < 5):
#          print "%12s -> %10s observed at angle %.2f on image %.1f" % (hkl, map_hkl_to_symmhkl[hkl], r['phi'],
#                                                                       expt.scan.get_image_index_from_angle(r['phi'] - (2 * 3.1415926535), deg=False))
      symmetry_mapped_hkls = set([map_hkl_to_symmhkl[hkl] for hkl in hkls])
      print "%5d unique reflections fall on detector during sweep (including symmetry relations)" % len(symmetry_mapped_hkls)
      completeness = 100 * len(symmetry_mapped_hkls) / completeness_limit
      multiplicity = len(hkls) / len(symmetry_mapped_hkls)
      score = evaluation_function(hkls)
      print "Estimated sweep completeness: %5.1f %%   sweep multiplicity : %.1f   sweep score: %.1f" %\
            (completeness, multiplicity, score)
      return {'completeness': completeness, 'multiplicity': multiplicity, 'score': score}

#    # count the number of observations per reflection to construct an evaluation function
#    seen_hkl_multiplicity = {x: 0 for x in possible_hkl}
#    for hkl in detectable_rays['miller_index']:
#      seen_hkl_multiplicity[hkl] += 1

    expt = experiments[0]
    spacegroup = expt.crystal.get_space_group()
    unit_cell = expt.crystal.get_unit_cell()

    possible_hkl = self.list_possible_reflections(spacegroup, unit_cell, dmin, dmax)
    possible_hkl_flex = flex.miller_index(possible_hkl)
    hkl_to_id = { hkl: id for (id, hkl) in enumerate(possible_hkl) }
    id_to_hkl = [ hkl for (id, hkl) in enumerate(possible_hkl) ]

    # find mapping of reciprocal space onto reciprocal asymmetric unit and its inverse
    from cctbx.miller import map_to_asu
    asu_hkl_flex = flex.miller_index(possible_hkl)
    map_to_asu(spacegroup.type(), False, asu_hkl_flex)
    # TODO: Treat anomalous signal?
    map_hkl_to_symmhkl = {r: rs for (r, rs) in zip(possible_hkl, list(asu_hkl_flex))}
    map_symmhkl_to_hkl = {}
    for k, v in map_hkl_to_symmhkl.iteritems():
      map_symmhkl_to_hkl[v] = map_symmhkl_to_hkl.get(v, [])
      map_symmhkl_to_hkl[v].append(k)

    unique_asu_indices = set(asu_hkl_flex)
    completeness_limit = len(unique_asu_indices)

    # scoring function: hkl_id -> shared_bin_id
    scoring_unique_reflection = { hkl: unique_id for (unique_id, hkl) in enumerate(possible_hkl) }
    scoring_unique_reflection = { hkl_to_id[h]: id for (h, id) in scoring_unique_reflection.iteritems() }
    scoring_symmetry_related_reflection = { hkl_to_id[hkl]: hkl_to_id[asu] for (hkl, asu) in map_hkl_to_symmhkl.iteritems() }

    scorings = { 'unique':   (scoring_unique_reflection,           [0] * len(possible_hkl)),
                 'symmetry': (scoring_symmetry_related_reflection, [0] * len(possible_hkl)),
               }

    print "%5d unique reflections possible in %s (ignoring anomalous signal)" %\
          (completeness_limit, spacegroup.type().lookup_symbol())
    print

    # Determine the detectable reflections given some data collection sweep
    s0 = expt.beam.get_s0()
    rotation_axis = expt.goniometer.get_rotation_axis()

    # Obtain the A matrix with all diffractometer angles set to 0.
    crystal_R = matrix.sqr(expt.goniometer.get_fixed_rotation())
    crystal_A = crystal_R.inverse() * expt.crystal.get_A()

    def run_strategies(strategies):
      # starting from scratch
      observed_hkls = {}

      total_score = 0.0
      degrees = 0

      from itertools import izip, count
      for (run, strategy) in izip(count(1), strategies):
        print
        print "Sweep %d:" % run
        evaluation_function = evaluation_function_factory(possible_hkl, observed_hkls, map_hkl_to_symmhkl, map_symmhkl_to_hkl)

        goniometer = goniometer_factory.make_kappa_goniometer(alpha=54.7,
                                                              kappa=strategy['kappa'],
                                                              phi=strategy['phi'],
                                                              omega=strategy['omega'],
                                                              direction="-y",
                                                              scan_axis=strategy['scan_axis'])

        degrees += strategy['scan']

        from math import radians
        if strategy['scan_axis'] == 'omega':
          oscillation_range = (radians(strategy['omega']), radians(strategy['omega'] + strategy['scan']))
        else:
          oscillation_range = (radians(strategy['phi']), radians(strategy['phi'] + strategy['scan']))

        # repeat sweep: expt.goniometer
        proposed_hkls = calculate_observations(expt.detector, goniometer, oscillation_range)

        strategy_results = evaluate_concrete_strategy(proposed_hkls, evaluation_function)
        combined_observations = self.add_seen_multiplicity(observed_hkls, proposed_hkls)

        num_observed_hkls      = len(set([map_hkl_to_symmhkl[hkl] for hkl in combined_observations.keys()]))
        count_hkl_observations = sum(combined_observations.values())
        completeness = 100 * num_observed_hkls / completeness_limit
        multiplicity = count_hkl_observations / num_observed_hkls
        total_score += strategy_results['score']
        print "Estimated total completeness: %5.1f %%   total multiplicity : %.1f   total score: %.1f" %\
              (completeness, multiplicity, total_score)

        # keep the repeat sweep and continue
        observed_hkls = combined_observations
      export_mtz(observed_hkls, expt, "mtzout.mtz")
      return {'completeness': completeness, 'multiplicity': multiplicity, 'score': total_score, 'degrees': degrees}

    return run_strategies


  def SimpleEvaluationFunctionFactory(self, possible_hkl, seen_hkl_multiplicity, map_hkl_to_symmhkl, map_symmhkl_to_hkl):
    # Asymptotic scoring.
    # First observation of a reflection: 1/1. Second: 1/2. Third: 1/3. etc.
    # Same scoring added again for the ASU reduced reflection observation
    #
    # This function returns a very basic evaluation function based on the number of observations per reflection
    # TODO: Evaluation function does not consider resolution shells of reflections

    reflection_count = {hkl: seen_hkl_multiplicity[hkl] if hkl in seen_hkl_multiplicity else 0 for hkl in possible_hkl}
    symm_reflection_count = {}
    for symm_hkl, symm_hkl_group in map_symmhkl_to_hkl.iteritems():
      symm_reflection_count[symm_hkl] = sum([reflection_count[hkl] for hkl in symm_hkl_group])

    def evaluation_function(hkls):
      _reflection_count = reflection_count.copy()
      _symm_reflection_count = symm_reflection_count.copy()

      score = 0.0
      for hkl in hkls:
        _reflection_count[hkl] += 1
        score += 1.0 / _reflection_count[hkl]
        _symm_reflection_count[map_hkl_to_symmhkl[hkl]] += 1
        score += 1.0 / _symm_reflection_count[map_hkl_to_symmhkl[hkl]]
      return score
    return evaluation_function

  def SimpleGeometricEvaluationFunctionFactory(self, possible_hkl, seen_hkl_multiplicity, map_hkl_to_symmhkl, map_symmhkl_to_hkl):
    # Geometric scoring.
    # First observation of a reflection: 0.5^0. Second: 0.5^1. Third: 0.5^2. etc.
    # Plus the same scoring scheme again for the ASU reduced reflection observation
    # Advantage: Order of observations does not matter for scoring purposes
    #
    # This function returns a very basic evaluation function based on the number of observations per reflection
    # TODO: Evaluation function does not consider resolution shells of reflections

    reflection_count = {hkl: seen_hkl_multiplicity[hkl] if hkl in seen_hkl_multiplicity else 0 for hkl in possible_hkl}
    symm_reflection_count = {}
    for symm_hkl, symm_hkl_group in map_symmhkl_to_hkl.iteritems():
      symm_reflection_count[symm_hkl] = sum([reflection_count[hkl] for hkl in symm_hkl_group])

    def evaluation_function(hkls):
      from math import pow
      _reflection_count = reflection_count.copy()
      _symm_reflection_count = symm_reflection_count.copy()

      score = 0.0
#      scores = []
      for hkl in hkls:
        score += pow(0.5, _reflection_count[hkl])
        score += pow(0.5, _symm_reflection_count[map_hkl_to_symmhkl[hkl]])
        _reflection_count[hkl] += 1
        _symm_reflection_count[map_hkl_to_symmhkl[hkl]] += 1
        #scores.append(score)

#      from collections import Counter
#      c = Counter(scores)
#      print c
#      return sum(scores)
      return score
    return evaluation_function

  def SimpleGeometricEvaluationFunctionFactory(self, possible_hkl, seen_hkl_multiplicity, map_hkl_to_symmhkl, map_symmhkl_to_hkl):
    # Geometric scoring.
    # First observation of a reflection: 0.5^0. Second: 0.5^1. Third: 0.5^2. etc.
    # Plus the same scoring scheme again for the ASU reduced reflection observation
    # Advantage: Order of observations does not matter for scoring purposes
    #
    # This function returns a very basic evaluation function based on the number of observations per reflection
    # TODO: Evaluation function does not consider resolution shells of reflections

    reflection_count = {hkl: seen_hkl_multiplicity[hkl] if hkl in seen_hkl_multiplicity else 0 for hkl in possible_hkl}
    symm_reflection_count = {}
    for symm_hkl, symm_hkl_group in map_symmhkl_to_hkl.iteritems():
      symm_reflection_count[symm_hkl] = sum([reflection_count[hkl] for hkl in symm_hkl_group])

    def evaluation_function(hkls):
      from math import pow
      _reflection_count = reflection_count.copy()
      _symm_reflection_count = symm_reflection_count.copy()

      score = 0.0
#      scores = []
      for hkl in hkls:
        score += pow(0.5, _reflection_count[hkl])
        score += pow(0.5, _symm_reflection_count[map_hkl_to_symmhkl[hkl]])
        _reflection_count[hkl] += 1
        _symm_reflection_count[map_hkl_to_symmhkl[hkl]] += 1
        #scores.append(score)

#      from collections import Counter
#      c = Counter(scores)
#      print c
#      return sum(scores)
      return score
    return evaluation_function

  def getFixedStrategylist(self):
    strategylist = []

    strategylist.append(
        {'name': "default sphere, K=75",
         'strategy': [
           {'kappa':  0, 'phi':   0, 'omega': -92, 'scan_axis': 'phi', 'scan': 360, '2theta': 30}

#           {'kappa': 75, 'phi': 160, 'omega': -102, 'scan_axis': 'omega', 'scan': 50, '2theta': 30} ,
#          {'kappa': 75, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
#          {'kappa': 75, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
#          {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    return strategylist

    strategylist.append(
        {'name': "20140302-runlist-sphere",
         'strategy': [
           {'kappa': 60, 'phi': 120, 'omega': -99, 'scan_axis': 'omega', 'scan': 134, '2theta': 25},
           {'kappa': 60, 'phi':   0, 'omega': -99, 'scan_axis': 'omega', 'scan': 134, '2theta': 25},
           {'kappa': 60, 'phi':-120, 'omega': -99, 'scan_axis': 'omega', 'scan': 134, '2theta': 25},
           {'kappa': 60, 'phi':-180, 'omega':-149, 'scan_axis': 'omega', 'scan': 134, '2theta':-25},
           {'kappa': 60, 'phi': -60, 'omega':-149, 'scan_axis': 'omega', 'scan': 133, '2theta':-25},
           {'kappa': 60, 'phi':  60, 'omega':-149, 'scan_axis': 'omega', 'scan': 133, '2theta':-25},
           {'kappa':  0, 'phi':  60, 'omega':-135, 'scan_axis': 'omega', 'scan': 180, '2theta': 25},
           {'kappa':  0, 'phi': -60, 'omega':-135, 'scan_axis': 'omega', 'scan': 180, '2theta': 25},
           {'kappa':  0, 'phi': -60, 'omega':  45, 'scan_axis': 'phi'  , 'scan': 180, '2theta': 25}
         ]
        }
    )

    strategylist.append(
        {'name': "cpRunlistAugust2013",
         'strategy': [
           {'kappa':-141, 'phi':  64, 'omega': -61, 'scan_axis': 'omega', 'scan': 136, '2theta':119},
           {'kappa':-107, 'phi': 149, 'omega': -61, 'scan_axis': 'omega', 'scan': 136, '2theta':119},
           {'kappa':-105, 'phi':-101, 'omega': -61, 'scan_axis': 'omega', 'scan': 114, '2theta':119},
           {'kappa': 109, 'phi': 117, 'omega':   1, 'scan_axis': 'omega', 'scan': 128, '2theta':119},
           {'kappa': 128, 'phi':-176, 'omega':  72, 'scan_axis': 'omega', 'scan':  60, '2theta':119},
           {'kappa': 144, 'phi':  20, 'omega':  24, 'scan_axis': 'omega', 'scan': 105, '2theta':119},
           {'kappa': 153, 'phi':-120, 'omega': -11, 'scan_axis': 'omega', 'scan': 140, '2theta':119}
         ]
        }
    )

#    return strategylist

    strategylist.append(
        {'name': "default sphere, K=65",
         'strategy': [
           {'kappa': 65, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 65, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 65, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "default sphere, K=75",
         'strategy': [
           {'kappa': 75, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 75, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 75, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "default sphere, K=90",
         'strategy': [
           {'kappa': 90, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 90, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 90, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "default sphere, K=60",
         'strategy': [
           {'kappa': 60, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 60, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 60, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "default sphere, K=45",
         'strategy': [
           {'kappa': 45, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 45, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 45, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa':  0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "Naive Phi",
         'strategy': [
           {'kappa': 0, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "140505-gascell",
         'strategy': [
           {'kappa': 0, 'phi': 160, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi':  40, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi': -80, 'omega': -92, 'scan_axis': 'omega', 'scan': 132, '2theta': 30},
           {'kappa': 0, 'phi': -80, 'omega':-160, 'scan_axis': 'omega', 'scan': 180, '2theta': 30}
         ]
        }
    )

    strategylist.append(
        {'name': "150615-standard",
         'strategy': [
           {'kappa': 0, 'phi':   0, 'omega': -92, 'scan_axis': 'phi',   'scan': 180, '2theta': 30},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'omega', 'scan': 180, '2theta': 30},
           {'kappa': 0, 'phi': 120, 'omega': -92, 'scan_axis': 'omega', 'scan': 180, '2theta': 30},
           {'kappa': 0, 'phi':  60, 'omega': -92, 'scan_axis': 'omega', 'scan': 180, '2theta': 30},
           {'kappa': 0, 'phi': 180, 'omega':   0, 'scan_axis': 'phi',   'scan': 180, '2theta': 30},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'omega', 'scan': 360, '2theta':-30},
           {'kappa': 0, 'phi': 180, 'omega':   0, 'scan_axis': 'phi',   'scan': 180, '2theta':-30},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'omega', 'scan': 360, '2theta': 60},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'omega', 'scan': 360, '2theta': 90},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'phi', 'scan': 180, '2theta': 60},
           {'kappa': 0, 'phi': 180, 'omega': -92, 'scan_axis': 'phi', 'scan': 180, '2theta': 90},
         ]
        }
    )

    return strategylist

  def run(self):
    """lurrr"""
    from dials.util.options import flatten_experiments

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True, args=None)
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) == 0:
      print "No experiments given"
      return
    elif len(experiments.imagesets()) > 1 or len(experiments.detectors()) > 1:
      raise RuntimeError('experiment list contains > 1 imageset or detector')

    params.prediction.dmin = 0.67
    params.prediction.dmax = 23
    if (params.prediction.dmin is None) or (params.prediction.dmax is None):
      integrated = pickle.load(open('/home/wra62962/local/testsets/Ni_dppe_NO2_2.dials.sweep3/integrated.pickle', 'rb'))
      #integrated = pickle.load(open('/home/wra62962/local/testsets/TestsetMX.300.dials/integrated.pickle', 'rb'))
      if params.prediction.dmin is None:
        params.prediction.dmin = min(integrated["d"])
      if params.prediction.dmax is None:
        params.prediction.dmax = max(integrated["d"])

    print "Considering reflections with resolution %.2f - %.2f Ang" % (params.prediction.dmin, params.prediction.dmax)

    evaluator = self.StrategyEvaluator(experiments,
                                       self.SimpleGeometricEvaluationFunctionFactory,
                                       #self.SimpleEvaluationFunctionFactory,
                                       params.prediction.dmin,
                                       params.prediction.dmax)

    strategylist = self.getFixedStrategylist()
    results = []

    for strategy in strategylist:
      print
      print "%s strategy: " % strategy['name']
      quality = evaluator(strategy['strategy'])
      quality['name'] = strategy['name']
      results.append(quality)

    print "%30s   Comp   Mul  Score  Sweep  Sc/deg" % "Strategy"
    for r in sorted(results, key=(lambda x: x['score']), reverse=True):
      print "%30s: %5.1f  %4.1f  %5d  %4d%s  %5.1f" %\
            (r['name'], r['completeness'], r['multiplicity'], r['score'], r['degrees'], " deg",
             r['score'] / r['degrees'])

if __name__ == '__main__':
  script = Script()
  script.run()
  sys.exit(0)

  from dxtbx.model.detector import detector_factory
  sensor = 'PAD'
  distance = 85
  beam_centre = (0, 0)
  pixel_size = (0.172, 0.172)
  image_size = (487, 619)

  two_theta_angle = 30
  d = detector_factory.two_theta(sensor, distance, beam_centre, '+x', '-y', '-x', two_theta_angle, pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [], px_mm = None)
  print d

  two_theta_angle = 10
  d = detector_factory.two_theta(sensor, distance, beam_centre, '+x', '-y', '-x', two_theta_angle, pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [], px_mm = None)
  print d

  for p in d:
    print p
    print dir(p)

    from scitbx import matrix
    # Obtain and normalize axes
    fast_axis = matrix.col(p.get_fast_axis()).normalize()
    slow_axis = matrix.col(p.get_slow_axis()).normalize()
    normal = fast_axis.cross(slow_axis)

    # Get the local transformation matrix
    tl = matrix.sqr(
     fast_axis.elems       + (0,) +
     slow_axis.elems       + (0,) +
     normal.elems          + (0,) +
     tuple(p.get_origin()) + (1,)).transpose()

    # Obtain 2theta=0 orientation
    two_theta_axis = matrix.col((-1.0, 0, 0))
    rinv = two_theta_axis.axis_and_angle_as_r3_rotation_matrix(-10, deg=True).transpose()
    rinv = matrix.sqr(
      rinv[0:3] + (0,) +
      rinv[3:6] + (0,) +
      rinv[6:9] + (0,) +
      (0, 0, 0, 1)).transpose()
    print "Undoing 2theta rotation:"
    print rinv
    print

    print "2theta=0 orientation matrix:"
    ztl = rinv * tl
    print ztl
    print

    # Obtain new 2theta orientation
    r = two_theta_axis.axis_and_angle_as_r3_rotation_matrix(30, deg=True).transpose()
    r = matrix.sqr(
      r[0:3] + (0,) +
      r[3:6] + (0,) +
      r[6:9] + (0,) +
      (0, 0, 0, 1)).transpose()

    tgt = (r * ztl).transpose()
    print "New 2theta orientation:"
    print tgt
    print

    p.set_frame(tgt[0:3], tgt[4:7], tgt[12:15])
    print p

  print "Theoretically equivalent to:"
  two_theta_angle = 30
  d = detector_factory.two_theta(sensor, distance, beam_centre, '+x', '-y', '-x', two_theta_angle, pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [], px_mm = None)
  print d
