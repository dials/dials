from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.array_family import flex

help_message = '''
'''

phil_scope = iotbx.phil.parse("""
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    check_format=False,
    epilog=help_message)

  params, options, args = parser.parse_args(show_diff_phil=True,
                                            return_unhandled=True)

  assert len(args) == 2
  from iotbx.reflection_file_reader import any_reflection_file

  xyz = []
  intensities = []
  lp_corrections = []

  for f in args:
    xdet = None
    ydet = None
    rot = None
    i_sigi = None
    lp = None
    arrays = any_reflection_file(f).as_miller_arrays(merge_equivalents=False)
    for ma in arrays:
      print ma.info().labels
      if ma.info().labels[0] == 'XDET':
        xdet = ma
      elif ma.info().labels[0] == 'YDET':
        ydet = ma
      elif ma.info().labels[0] == 'ROT':
        rot = ma
      elif ma.info().labels == ['I', 'SIGI']:
        i_sigi = ma
      elif ma.info().labels[0] == 'LP':
        lp = ma

    assert [xdet, ydet, rot, i_sigi, lp].count(None) == 0

    xyz.append(flex.vec3_double(xdet.data(), ydet.data(), rot.data()))
    intensities.append(i_sigi)
    lp_corrections.append(lp)

  xyz1, xyz2 = xyz
  xyz2 += (1e-3,1e-3,1e-3)
  intensities1, intensities2 = intensities
  lp1, lp2 = lp_corrections

  # Do the nn match
  from annlib_ext import AnnAdaptor as ann_adaptor
  ann = ann_adaptor(xyz1.as_double().as_1d(), 3)
  ann.query(xyz2.as_double().as_1d())

  distances = flex.sqrt(ann.distances)
  matches = distances < 2 #pixels
  index1 = flex.size_t(list(ann.nn.select(matches)))
  index2 = flex.size_t(list(matches.iselection()))

  intensities1 = intensities1.select(index1)
  intensities2 = intensities2.select(index2)
  isigi1 = intensities1.data()/intensities1.sigmas()
  isigi2 = intensities2.data()/intensities2.sigmas()
  lp1 = lp1.select(index1)
  lp2 = lp2.select(index2)
  ##differences = intensities1.data() - intensities2.data()
  ##sums = intensities1.data() + intensities2.data()
  #differences = isigi1 - isigi2
  #sums = isigi1 + isigi2
  #assert sums.all_ne(0)
  #dos = differences/sums

  #mean_dos = []
  #binner = intensities1.setup_binner_d_star_sq_step(d_star_sq_step=0.01)
  #d_spacings = intensities1.d_spacings().data()
  #for i in range(binner.n_bins_used()):
    #d_max, d_min = binner.bin_d_range(i+1)
    #bin_sel = (d_spacings > d_min) & (d_spacings <= d_max)
    #mean_dos.append(flex.mean(dos.select(bin_sel)))

  # set backend before importing pyplot
  import matplotlib
  matplotlib.use('Agg')

  from matplotlib import pyplot
  pyplot.scatter(intensities1.data(), intensities2.data(), marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_intensities.png')
  pyplot.clf()

  pyplot.scatter(intensities1.sigmas(), intensities2.sigmas(), marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_sigmas.png')
  pyplot.clf()

  pyplot.scatter(
    flex.pow2(intensities1.sigmas()), flex.pow2(intensities2.sigmas()),
    marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_variances.png')
  pyplot.clf()

  pyplot.scatter(isigi1, isigi2, marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_i_sig_i.png')
  pyplot.clf()

  pyplot.scatter(lp1.data(), lp2.data(), marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m))
  pyplot.savefig('scatter_LP.png')
  pyplot.clf()

  #from cctbx import uctbx
  #pyplot.scatter(uctbx.d_star_sq_as_d(binner.bin_centers(2)), mean_dos)
  #pyplot.savefig('mean_dos.png')
  #pyplot.clf()


  return



def mean_difference_over_sum(self, other, binner=None):
  """
  """
  if (binner is not None):

    return min(self.indices().size() / max(1, complete_set.indices().size()),
               1.0) * multiplier
  assert self.binner() is not None
  data = []
  for n_given,n_complete in zip(self.binner().counts_given(),
                                self.binner().counts_complete()):
    if (n_complete == 0): data.append(return_fail)
    else: data.append(multiplier*n_given/n_complete)
  return binned_data(binner=self.binner(), data=data, data_fmt="%5.3f")



if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
