from __future__ import division
import copy
import math
from rstbx.symmetry.subgroup import MetricSubgroup
from scitbx.matrix import col
from libtbx import easy_mp


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
        'unit_cell':uc.parameters()
        }

    return result

  def labelit_printout(self,out=None):
    from libtbx import table_utils
    if out is None:
      import sys
      out = sys.stdout

    table_data = [["Solution","Metric fit","rmsd","#spots",
                   "crystal_system","unit_cell","volume"]]
    for item in self:
      uc = item.refined_crystal.get_unit_cell()
      P = uc.parameters()
      table_data.append(['%6d'%item.setting_number,
                         "%(max_angular_difference)6.4f dg"%item,
                         "%5.3f"%item.rmsd,
                         "%d"%item.Nmatches,
                         "%(system)s %(bravais)s"%item,
                         "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f"%P,
                         "%.0f"%uc.volume()])

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
    space_group = sgtbx.space_group_info(
      number=bravais_lattice_to_lowest_symmetry_spacegroup_number[bravais]).group()
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
  refinery, refined, outliers = refine(
    params, used_reflections, experiments, verbosity=refiner_verbosity)

  dall = refinery.rmsds()
  dx = dall[0]; dy = dall[1]
  subgroup.rmsd = math.sqrt(dx*dx + dy*dy)
  subgroup.Nmatches = len(refinery.get_matches())
  subgroup.scan = refinery.get_scan()
  subgroup.goniometer = refinery.get_goniometer()
  subgroup.beam = refinery.get_beam()
  subgroup.detector = refinery.get_detector()
  subgroup.refined_crystal = refinery.get_crystal()
  return subgroup
