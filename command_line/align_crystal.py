# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
import copy
from scitbx import matrix
import iotbx.phil
from cctbx import sgtbx

help_message = '''

'''

phil_scope= iotbx.phil.parse("""
space_group = None
  .type = space_group
align {
  mode = *main cusp
    .type = choice
  crystal {
    vector = None
      .type = str
      .multiple = True
    frame = *reciprocal direct
      .type = choice
  }
}
""")


a_star = matrix.col((1.,0.,0.))
b_star = matrix.col((0.,1.,0.))
c_star = matrix.col((0.,0.,1.))

a = matrix.col((1.,0.,0.))
b = matrix.col((0.,1.,0.))
c = matrix.col((0.,0.,1.))


class align_crystal(object):

  vector_names = {
    a.elems: 'a',
    b.elems: 'b',
    c.elems: 'c',
  }

  def __init__(self, experiment, vectors, frame='reciprocal', mode='main'):
    self.experiment = experiment
    self.vectors = vectors
    self.frame = frame
    self.mode = mode

    gonio = experiment.goniometer
    scan = experiment.scan

    self.s0 = matrix.col(self.experiment.beam.get_s0())
    self.rotation_axis = matrix.col(gonio.get_rotation_axis())

    axes = gonio.get_axes()
    assert len(axes) == 3
    e1, e2, e3 = (matrix.col(e) for e in axes)

    fixed_rotation = matrix.sqr(gonio.get_fixed_rotation())
    setting_rotation = matrix.sqr(gonio.get_setting_rotation())
    rotation_axis = matrix.col(gonio.get_rotation_axis())
    rotation_matrix = rotation_axis.axis_and_angle_as_r3_rotation_matrix(
      experiment.scan.get_oscillation()[0], deg=True)

    from dials.algorithms.refinement import rotation_decomposition

    results = {}

    for (v1_, v2_) in self.vectors:
      results[(v1_, v2_)] = {}
      crystal = copy.deepcopy(self.experiment.crystal)
      for smx in list(crystal.get_space_group().smx())[:]:
        cb_op = sgtbx.change_of_basis_op(smx)
        crystal = crystal.change_basis(cb_op)

        # Goniometer datum setting [D] at which the orientation was determined
        D = (setting_rotation * rotation_matrix * fixed_rotation).inverse()

        # The setting matrix [U] will vary with the datum setting according to
        # [U] = [D] [U0]
        U = crystal.get_U()

        # XXX In DIALS recorded U is equivalent to U0 - D is applied to U inside
        # prediction
        U0 = U

        B = crystal.get_B()

        if self.frame == 'direct':
          B = B.inverse().transpose()

        v1_0 = U0 * B * v1_
        v2_0 = U0 * B * v2_

        #c  (b) The laboratory frame vectors l1 & l2 are normally specified with the
        #c MODE command: MODE MAIN (the default) sets l1 (along which v1 will be
        #c placed) along the principle goniostat axis e1 (Omega), and l2 along
        #c the beam s0. This allows rotation for instance around a principle axis.
        #c The other mode is MODE CUSP, which puts l1 (v1) perpendicular to the
        #c beam (s0) and the e1 (Omega) axis, and l2 (v2) in the plane containing
        #c l1 & e1 (ie l1 = e1 x s0, l2 = e1).

        if self.mode == 'cusp':
          l1 = self.rotation_axis.cross(s0)
          l2 = self.rotation_axis
        else:
          l1 = self.rotation_axis.normalize()
          l3 = l1.cross(self.s0).normalize()
          l2 = l1.cross(l3)

        from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
        R = align_reference_frame(v1_0, l1, v2_0, l2)

        solutions = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
          R, e1, e2, e3, return_both_solutions=True, deg=True)

        if solutions is None:
          continue

        results[(v1_, v2_)][smx] = solutions

    self.all_solutions = results

    self.unique_solutions = {}
    for (v1, v2), result in results.iteritems():
      for solutions in result.itervalues():
        for solution in solutions:
          k = tuple(round(a, 2) for a in solution[1:])
          self.unique_solutions.setdefault(k, set())
          self.unique_solutions[k].add((v1, v2))

  def show(self):
    from libtbx import table_utils
    self.info()

    rows = []
    names = self.experiment.goniometer.get_names()
    rows.append(['v1', 'v2', names[0], names[1]])

    def vector_as_str(v):
      v = v.elems
      if v in self.vector_names:
        vstr = self.vector_names[v]
        if self.frame == 'reciprocal':
          vstr += '*'
      else:
        vstr = str(v)
      return vstr

    for angles, solutions in self.unique_solutions.iteritems():
        for (v1, v2) in solutions:
          rows.append(
            (vector_as_str(v1), vector_as_str(v2),
             '% 7.2f' %angles[0], '% 7.2f' %angles[1],
             ))
    print 'Independent solutions:'
    print table_utils.format(rows=rows, has_header=True)

  def info(self):
    from libtbx import table_utils

    U = self.experiment.crystal.get_U()
    B = self.experiment.crystal.get_B()

    a_star_ = U * B * a_star
    b_star_ = U * B * b_star
    c_star_ = U * B * c_star

    Binvt = B.inverse().transpose()

    a_ = U * Binvt * a
    b_ = U * Binvt * b
    c_ = U * Binvt * c

    def smallest_angle(angle):
      return min(angle, 180-angle)

    rows = [['Experimental axis', 'a*', 'b*', 'c*']]
    rows.append(['Rotation'] + [
      '%.3f' %smallest_angle(axis.angle(self.rotation_axis, deg=True))
      for axis in (a_star_, b_star_, c_star_)])
    rows.append(['Beam'] + [
      '%.3f' %smallest_angle(axis.angle(self.s0, deg=True))
      for axis in (a_star_, b_star_, c_star_)])
    print 'Angles between reciprocal cell axes and principal experimental axes:'
    print table_utils.format(rows=rows, has_header=True)
    print

    rows = [['Experimental axis', 'a', 'b', 'c']]
    rows.append(['Rotation'] + [
      '%.3f' %smallest_angle(axis.angle(self.rotation_axis, deg=True))
      for axis in (a_, b_, c_)])
    rows.append(['Beam'] + [
      '%.3f' %smallest_angle(axis.angle(self.s0, deg=True))
      for axis in (a_, b_, c_)])
    print 'Angles between unit cell axes and principal experimental axes:'
    print table_utils.format(rows=rows, has_header=True)
    print


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)

  if len(experiments) == 0:
    parser.print_help()
    exit(0)

  imagesets = experiments.imagesets()

  expt = experiments[0]

  frame = None

  if len(params.align.crystal.vector):
    assert len(params.align.crystal.vector) %2 == 0
    vectors = []

    name_to_vectors = {
      'a': (a, 'direct'),
      'b': (b, 'direct'),
      'c': (c, 'direct'),
      'a*': (a_star, 'reciprocal'),
      'b*': (b_star, 'reciprocal'),
      'c*': (c_star, 'reciprocal'),
    }

    for v in params.align.crystal.vector:
      v = v.strip()
      if v in name_to_vectors:
        v, frame_ = name_to_vectors[v]
        assert frame is None or frame == frame_
        frame = frame_
      else:
        v = v.replace(',', ' ').strip().split()
        assert len(v) == 3
        v = matrix.col([float(v_) for v_ in v])
        if frame is None:
          frame = params.align.crystal.frame

      vectors.append(v)
    vectors = [(vectors[2*i], vectors[2*i+1])
               for i in range(len(vectors)//2)]
  else:
    frame = 'reciprocal'
    vectors = ((a_star, b_star), # a*, b*
               (a_star, c_star), # a*, c*
               (b_star, a_star), # b*, a*
               (b_star, c_star), # b*, c*
               (c_star, a_star), # c*, a*
               (c_star, b_star), # c*, b*
              )


  result = align_crystal(expt, vectors, frame=frame, mode=params.align.mode)
  result.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
