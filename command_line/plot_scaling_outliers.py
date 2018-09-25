"""
Plot the outliers determined during scaling from a scaled.pickle
"""

from __future__ import absolute_import, division, print_function
import sys
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections
from libtbx import phil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



phil_scope = phil.parse('''
  output {
    plot_out = "outliers.png"
      .type = str
      .help = "Option to set filename for output plot."
  }
''')

def main(argv):
  '''the plotting script'''

  optionparser = OptionParser(usage=None, read_experiments=True,
    read_reflections=True, phil=phil_scope,
    check_format=False)
  params, _ = optionparser.parse_args(argv, show_diff_phil=False)


  if not params.input.reflections:
    optionparser.print_help()
    return

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    print('The following parameters have been modified:\n')
    print(diff_phil)

  reflections = flatten_reflections(params.input.reflections)[0]
  if len(set(reflections['id'])) == 1:
    plot_outliers(params, reflections, "outliers.png")
  else:
    ids = set(reflections['id'])
    for i, id_val in enumerate(ids):
      table = reflections.select(reflections['id'] == id_val)
      plot_outliers(params, table, filename="outliers_"+str(i+1)+".png")

def plot_outliers(params, reflections, filename=None):
  '''plots positions of outliers'''
  outliers = reflections.select(reflections.get_flags(
    reflections.flags.outlier_in_scaling))
  plt.figure(figsize=(11, 6))
  gs = gridspec.GridSpec(1, 2)
  ax1 = plt.subplot(gs[0, 0])
  ax2 = plt.subplot(gs[0, 1])
  x, y, z = outliers['xyzobs.px.value'].parts()
  n_bins = 100
  ax1.scatter(x, y, s=0.5)
  ax1.set_ylabel('detector y-position')
  ax1.set_xlabel('detector x-position')
  ax2.hist(z, bins=n_bins, label="n_bins = %s" % n_bins)
  ax2.set_ylabel('number of outliers')
  ax2.set_xlabel('image number')
  ax2.legend()

  if filename:
    print("Saving plot to %s" % filename)
    plt.savefig(filename)
  else:
    print("Saving plot to %s" % params.output.plot_out)
    plt.savefig(params.output.plot_out)

if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
