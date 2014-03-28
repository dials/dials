#!/usr/bin/env python
#
# import.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import matplotlib

# Offline backend
matplotlib.use("Agg")

from matplotlib import pylab

def ensure_directory(path):
  ''' Make the directory if not already there. '''
  from os import makedirs
  import errno
  try:
    makedirs(path)
  except OSError as e:
    if e.errno != errno.EEXIST:
      raise

def ensure_required(rlist, required):
  ''' Check which keys aren't present. '''
  not_present = []
  for k in required:
    if k not in rlist:
      not_present.append(k)
  if len(not_present) != 0:
    print " Skipping: following required fields not present:"
    for k in not_present:
      print "  %s" % k
    return False
  return True


class ReferenceProfileAnalyser(object):
  ''' Analyse the reference profiles. '''

  def __init__(self, directory):
    ''' Set up the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "reference")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.raw.value",
      "intensity.raw.variance",
      "xyzcal.px",
      "profile.correlation"
    ]

  def __call__(self, rlist):
    ''' Analyse the reference profiles. '''

    # Check we have the required fields
    print "Analysing reference profiles"
    if not ensure_required(rlist, self.required):
      return

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    print " Analysing reference profile distribution vs x/y"
    self.reference_xy(rlist)

    print " Analysing reference profile distribution vs z"
    self.reference_z(rlist)

    def correlations(filename, rlist):
      ''' Call for reference spots and all reflections. '''

      print " Analysing reflection profile correlations"
      self.reflection_corr_hist(rlist, filename)

      print " Analysing reflection profile correlations vs x/y"
      self.reflection_corr_vs_xy(rlist, filename)

      print " Analysing reflection profile correlations vs z"
      self.reflection_corr_vs_z(rlist, filename)

      print " Analysing reflection profile correlations vs I/Sigma"
      self.reflection_corr_vs_ios(rlist, filename)

    mask = rlist.get_flags(rlist.flags.reference_spot)
    correlations("reference",  rlist.select(mask))
    correlations("reflection", rlist)

  def reference_xy(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from dials.array_family import flex
    from os.path import join
    mask = rlist.get_flags(rlist.flags.reference_spot)
    rlist = rlist.select(mask)
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reference profiles binned in X/Y")
    cax = pylab.hexbin(x, y, gridsize=20)
    pylab.colorbar(cax)
    pylab.savefig(join(self.directory, "reference_xy.png"))
    pylab.clf()

  def reference_z(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations vs Z")
    cax = pylab.hist(z, bins=20)
    pylab.savefig(join(self.directory, "reference_z.png"))
    pylab.clf()

  def reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    pylab.title("Reflection correlations histogram")
    pylab.hist(corr, bins=20)
    pylab.savefig(join(self.directory, "%s_corr_hist" % filename))
    pylab.clf()

  def reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations binned in X/Y")
    cax = pylab.hexbin(x, y, C=corr, gridsize=20, vmin=0.0, vmax=1.0)
    pylab.colorbar(cax)
    pylab.savefig(join(self.directory, "%s_corr_vs_xy.png" % filename))
    pylab.clf()

  def reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations vs Z")
    cax = pylab.hexbin(z, corr, gridsize=20)
    pylab.colorbar(cax)
    pylab.savefig(join(self.directory, "%s_corr_vs_z.png" % filename))
    pylab.clf()

  def reflection_corr_vs_ios(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > 0.1
    I_over_S = I_over_S.select(mask)
    corr = corr.select(mask)
    pylab.title("Reflection correlations vs Log I/Sigma")
    cax = pylab.hexbin(flex.log(I_over_S), corr, gridsize=19)
    pylab.colorbar(cax)
    pylab.savefig(join(self.directory, "%s_corr_vs_ios.png" % filename))
    pylab.clf()


class Analyser(object):
  ''' Helper class to do all the analysis. '''

  def __init__(self, directory):
    ''' Setup the analysers. '''
    from os.path import join
    directory = join(directory, "analysis")
    self.analysers = [
      ReferenceProfileAnalyser(directory),
    ]

  def __call__(self, rlist):
    ''' Do all the analysis. '''
    for analyse in self.analysers:
      analyse(rlist)


if __name__ == '__main__':

  from dials.array_family import flex
  from optparse import OptionParser
  from dials.util.command_line import Command
  usage = "usage: %prog [options] reflections.pickle"
  parser = OptionParser(usage)

  # The directory to store output
  parser.add_option(
    "-d", "--dir",
    dest = "dir",
    type = "string", default = ".",
    help = "The directory to store results")

  # Parse the command line arguments
  options, args = parser.parse_args()

  # Shoe the help
  if len(args) != 1:
    parser.show_help()
    exit(0)

  # Read the reflections
  Command.start("Loading reflections")
  rlist = flex.reflection_table.from_pickle(args[0])
  Command.end("Loaded %d reflections" % len(rlist))

  # Analyse the reflections
  analyse = Analyser(options.dir)
  analyse(rlist)
