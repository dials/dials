from __future__ import division
import copy
import math
from libtbx import easy_mp
from cctbx import crystal, sgtbx
from scitbx.matrix import col
from scitbx.array_family import flex
from rstbx.symmetry.subgroup import MetricSubgroup



def dials_crystal_from_orientation(crystal_orientation,space_group):
  dm = crystal_orientation.direct_matrix()
  AA = col((dm[0],dm[1],dm[2]))
  BB = col((dm[3],dm[4],dm[5]))
  CC = col((dm[6],dm[7],dm[8]))

  from dxtbx.model.crystal import crystal_model

  cryst = crystal_model(real_space_a=AA, real_space_b=BB, real_space_c=CC,
                        space_group=space_group)
  return cryst


class bravais_setting(MetricSubgroup): # inherits from dictionary
  def __init__(self,other):
    self.update(other)


class refined_settings_list(list):

  def supergroup(self):
    return self[0]

  def triclinic(self):
    return self[-1]

  def as_dict(self):
    result = { }

    for item in self:
      uc = item.refined_crystal.get_unit_cell()
      result[item.setting_number] = {
        'max_angular_difference':item['max_angular_difference'],
        'rmsd':item.rmsd,
        'nspots':item.Nmatches,
        'bravais':item['bravais'],
        'unit_cell':uc.parameters(),
        'cb_op':item['cb_op_inp_best'].as_abc(),
        'max_cc':item.max_cc,
        'min_cc':item.min_cc,
        }

    return result

  def labelit_printout(self,out=None):
    from libtbx import table_utils
    if out is None:
      import sys
      out = sys.stdout

    table_data = [["Solution","Metric fit","rmsd", "min/max cc", "#spots",
                   "lattice","unit_cell","volume", "cb_op"]]
    for item in self:
      uc = item.refined_crystal.get_unit_cell()
      P = uc.parameters()
      min_max_cc_str = "-/-"
      if item.min_cc is not None and item.max_cc is not None:
        min_max_cc_str = "%.3f/%.3f" %(item.min_cc, item.max_cc)
      table_data.append(['%6d'%item.setting_number,
                         "%(max_angular_difference)6.4f"%item,
                         "%5.3f"%item.rmsd,
                         min_max_cc_str,
                         "%d"%item.Nmatches,
                         "%(bravais)s"%item,
                         "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f"%P,
                         "%.0f"%uc.volume(),
                         "%s"%item['cb_op_inp_best'].as_abc()])

    print >> out, table_utils.format(
        table_data, has_header=1, justify='right', delim=' ')

# Mapping of Bravais lattice type to corresponding lowest possible symmetry
bravais_lattice_to_lowest_symmetry_spacegroup_number = {
  'aP':1, 'mP':3, 'mC':5, 'oP':16, 'oC':20, 'oF':22, 'oI':23, 'tP':75,
  'tI':79, 'hP':143, 'hR':146, 'cP':195, 'cF':196, 'cI':197
}

def refined_settings_factory_from_refined_triclinic(
  params, experiment, reflections, i_setting=None,
  lepage_max_delta=5.0, nproc=1, refiner_verbosity=0):

  crystal = experiment.crystal

  used_reflections = copy.deepcopy(reflections)
  UC = crystal.get_unit_cell()

  from rstbx.dps_core.lepage import iotbx_converter

  Lfat = refined_settings_list()
  for item in iotbx_converter(UC, lepage_max_delta):
    Lfat.append(bravais_setting(item))

  supergroup = Lfat.supergroup()
  triclinic = Lfat.triclinic()
  triclinic_miller = used_reflections['miller_index']

  # assert no transformation between indexing and bravais list
  assert str(triclinic['cb_op_inp_best'])=="a,b,c"

  Nset = len(Lfat)
  for j in xrange(Nset):  Lfat[j].setting_number = Nset-j

  from cctbx.crystal_orientation import crystal_orientation
  from cctbx import sgtbx
  from scitbx import matrix
  for j in xrange(Nset):
    cb_op = Lfat[j]['cb_op_inp_best'].c().as_double_array()[0:9]
    orient = crystal_orientation(crystal.get_A(),True)
    orient_best = orient.change_basis(matrix.sqr(cb_op).transpose())
    constrain_orient = orient_best.constrain(Lfat[j]['system'])
    bravais = Lfat[j]["bravais"]
    cb_op_best_ref = Lfat[j]['best_subsym'].change_of_basis_op_to_reference_setting()
    space_group = sgtbx.space_group_info(
      number=bravais_lattice_to_lowest_symmetry_spacegroup_number[bravais]).group()
    space_group = space_group.change_basis(cb_op_best_ref.inverse())
    bravais = str(bravais_types.bravais_lattice(group=space_group))
    Lfat[j]["bravais"] = bravais
    Lfat[j].unrefined_crystal = dials_crystal_from_orientation(
      constrain_orient, space_group)

  args = []
  for subgroup in Lfat:
    args.append((
      params, subgroup, used_reflections, experiment, refiner_verbosity))

  results = easy_mp.parallel_map(
    func=refine_subgroup,
    iterable=args,
    processes=nproc,
    method="multiprocessing",
    preserve_order=True,
    asynchronous=True,
    preserve_exception_message=True)

  for i, result in enumerate(results):
    Lfat[i] = result
  return Lfat

def refine_subgroup(args):
  assert len(args) == 5
  from dials.command_line.check_indexing_symmetry \
       import get_symop_correlation_coefficients

  params, subgroup, used_reflections, experiment, refiner_verbosity = args

  used_reflections = copy.deepcopy(used_reflections)
  triclinic_miller = used_reflections['miller_index']
  cb_op = subgroup['cb_op_inp_best']
  higher_symmetry_miller = cb_op.apply(triclinic_miller)
  used_reflections['miller_index'] = higher_symmetry_miller
  experiment.crystal = copy.deepcopy(subgroup.unrefined_crystal)

  from dials.algorithms.indexing.refinement import refine
  from dxtbx.model.experiment.experiment_list import ExperimentList
  experiments = ExperimentList([experiment])
  subgroup.max_cc = None
  subgroup.min_cc = None
  try:
    refinery, refined, outliers = refine(
      params, used_reflections, experiments, verbosity=refiner_verbosity)
  except RuntimeError, e:
    if (str(e) == "scitbx Error: g0 - astry*astry -astrz*astrz <= 0." or
        str(e) == "scitbx Error: g1-bstrz*bstrz <= 0."):
      subgroup.refined_crystal = None
      subgroup.rmsd = None
      subgroup.Nmatches = None
      subgroup.scan = None
      subgroup.goniometer = None
      subgroup.beam = None
      subgroup.detector = None
    else: raise
  else:
    dall = refinery.rmsds()
    dx = dall[0]; dy = dall[1]
    subgroup.rmsd = math.sqrt(dx*dx + dy*dy)
    subgroup.Nmatches = len(refinery.get_matches())
    subgroup.scan = refinery.get_scan()
    subgroup.goniometer = refinery.get_goniometer()
    subgroup.beam = refinery.get_beam()
    subgroup.detector = refinery.get_detector()
    subgroup.refined_crystal = refinery.get_crystal()
    cs = crystal.symmetry(
        unit_cell=subgroup.refined_crystal.get_unit_cell(),
        space_group=subgroup.refined_crystal.get_space_group())
    from cctbx import miller
    ms = miller.set(cs, used_reflections['miller_index'])
    ms = ms.array(used_reflections['intensity.sum.value'] /
                flex.sqrt(used_reflections['intensity.sum.variance']))
    ccs, nrefs = get_symop_correlation_coefficients(ms)
    subgroup.correlation_coefficients = ccs
    subgroup.cc_nrefs = nrefs
    if len(ccs) > 1:
      subgroup.max_cc = flex.max(ccs[1:])
      subgroup.min_cc = flex.min(ccs[1:])
  return subgroup

from cctbx.sgtbx import subgroups
from cctbx.sgtbx import bravais_types
from cctbx.sgtbx import change_of_basis_op
find_max_delta = sgtbx.lattice_symmetry_find_max_delta

def metric_supergroup(group):
  return sgtbx.space_group_info(group=group).type(
    ).expand_addl_generators_of_euclidean_normalizer(True,True
    ).build_derived_acentric_group()

def find_matching_symmetry(unit_cell, target_space_group, max_delta=5):
  cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())
  target_bravais_t = bravais_types.bravais_lattice(
    group=target_space_group.info().reference_setting().group())
  best_subgroup = None
  best_angular_difference = 1e8

  # code based on cctbx/sgtbx/lattice_symmetry.py but optimised to only
  # look at subgroups with the correct bravais type

  input_symmetry = cs
  # Get cell reduction operator
  cb_op_inp_minimum = input_symmetry.change_of_basis_op_to_minimum_cell()

  # New symmetry object with changed basis
  minimum_symmetry = input_symmetry.change_basis(cb_op_inp_minimum)

  # Get highest symmetry compatible with lattice
  lattice_group = sgtbx.lattice_symmetry_group(
    minimum_symmetry.unit_cell(),
    max_delta=max_delta,
    enforce_max_delta_for_generated_two_folds=True)

  # Get list of sub-spacegroups
  subgrs = subgroups.subgroups(lattice_group.info()).groups_parent_setting()

  # Order sub-groups
  sort_values = flex.double()
  for group in subgrs:
    order_z = group.order_z()
    space_group_number = sgtbx.space_group_type(group, False).number()
    assert 1 <= space_group_number <= 230
    sort_values.append(order_z*1000+space_group_number)
  perm = flex.sort_permutation(sort_values, True)

  for i_subgr in perm:
    acentric_subgroup = subgrs[i_subgr]
    acentric_supergroup = metric_supergroup(acentric_subgroup)
    ## Add centre of inversion to acentric lattice symmetry
    #centric_group = sgtbx.space_group(acentric_subgroup)
    #centric_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    # Make symmetry object: unit-cell + space-group
    # The unit cell is potentially modified to be exactly compatible
    # with the space group symmetry.
    subsym = crystal.symmetry(
      unit_cell=minimum_symmetry.unit_cell(),
      space_group=acentric_subgroup,
      assert_is_compatible_unit_cell=False)
    #supersym = crystal.symmetry(
      #unit_cell=minimum_symmetry.unit_cell(),
      #space_group=acentric_supergroup,
      #assert_is_compatible_unit_cell=False)
    # Convert subgroup to reference setting
    cb_op_minimum_ref = subsym.space_group_info().type().cb_op()
    ref_subsym = subsym.change_basis(cb_op_minimum_ref)
    # Ignore unwanted groups
    bravais_t = bravais_types.bravais_lattice(
      group=ref_subsym.space_group())
    if bravais_t != target_bravais_t:
      continue

    # Choose best setting for monoclinic and orthorhombic systems
    cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell(
      best_monoclinic_beta=True)

    best_subsym = ref_subsym.change_basis(cb_op_best_cell)
    # Total basis transformation
    cb_op_best_cell = change_of_basis_op(str(cb_op_best_cell),stop_chars='',r_den=144,t_den=144)
    cb_op_minimum_ref=change_of_basis_op(str(cb_op_minimum_ref),stop_chars='',r_den=144,t_den=144)
    cb_op_inp_minimum=change_of_basis_op(str(cb_op_inp_minimum),stop_chars='',r_den=144,t_den=144)
    cb_op_inp_best = cb_op_best_cell * cb_op_minimum_ref * cb_op_inp_minimum
    # Use identity change-of-basis operator if possible
    if (best_subsym.unit_cell().is_similar_to(input_symmetry.unit_cell())):
      cb_op_corr = cb_op_inp_best.inverse()
      try:
        best_subsym_corr = best_subsym.change_basis(cb_op_corr)
      except RuntimeError, e:
        if (str(e).find("Unsuitable value for rational rotation matrix.") < 0):
          raise
      else:
        if (best_subsym_corr.space_group() == best_subsym.space_group()):
          cb_op_inp_best = cb_op_corr * cb_op_inp_best

    max_angular_difference = find_max_delta(
      reduced_cell=minimum_symmetry.unit_cell(),
      space_group=acentric_supergroup)

    if max_angular_difference < best_angular_difference:
      #best_subgroup = subgroup
      best_angular_difference = max_angular_difference
      best_subgroup = {'subsym':subsym,
                       #'supersym':supersym,
                       'ref_subsym':ref_subsym,
                       'best_subsym':best_subsym,
                       'cb_op_inp_best':cb_op_inp_best,
                       'max_angular_difference':max_angular_difference
                       }

  if best_subgroup is not None:
    return best_subgroup
