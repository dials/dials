#!/usr/bin/env python
#
# analyse_output.py
#
#  Copyright (C) 2014 Diamond Light Source
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


class CentroidAnalyser(object):
  ''' Analyse the reflection centroids. '''

  def __init__(self, directory):
    ''' Setup the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "centroid")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.raw.value",
      "intensity.raw.variance",
      "xyzcal.px",
      "xyzobs.px.value"
    ]

  def __call__(self, rlist):
    ''' Analyse the reflection centroids. '''

    # Check we have the required fields
    print "Analysing reflection centroids"
    if not ensure_required(rlist, self.required):
      return

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    # Look at differences in calculated/observed position
    print " Analysing centroid differences with I/Sigma > 10"
    self.centroid_diff_hist(rlist, 10)
    print " Analysing centroid differences in x/y with I/Sigma > 10"
    self.centroid_diff_xy(rlist, 10)
    print " Analysing centroid differences in z with I/Sigma > 10"
    self.centroid_diff_z(rlist, 10)

  def centroid_diff_hist(self, rlist, threshold):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    xd = xo - xc
    yd = yo - yc
    zd = zo - zc
    diff = flex.sqrt(xd*xd + yd*yd + zd*zd)
    pylab.title("Difference between observed and calculated")
    cax = pylab.hist(diff, bins=20)
    pylab.xlabel("Difference in position")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "centroid_diff_hist.png"))
    pylab.clf()

  def centroid_diff_xy(self, rlist, threshold):
    ''' Look at the centroid difference in x, y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    xd = xo - xc
    yd = yo - yc
    pylab.title("Difference between observed and calculated in X")
    cax = pylab.hexbin(xc, yc, C=xd, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Difference in x position")
    pylab.savefig(join(self.directory, "centroid_diff_x.png"))
    pylab.clf()
    pylab.title("Difference between observed and calculated in Y")
    cax = pylab.hexbin(xc, yc, C=yd, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Difference in y position")
    pylab.savefig(join(self.directory, "centroid_diff_y.png"))
    pylab.clf()

  def centroid_diff_z(self, rlist, threshold):
    ''' Look at the centroid difference in x, y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    mask = I_over_S > threshold
    rlist = rlist.select(mask)
    assert(len(rlist) > 0)
    xc, yc, zc = rlist['xyzcal.px'].parts()
    xo, yo, zo = rlist['xyzobs.px.value'].parts()
    zd = zo - zc
    pylab.title("Difference between observed and calculated in Z")
    cax = pylab.hexbin(zc, zd, gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Difference in z position")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# Reflections")
    pylab.savefig(join(self.directory, "centroid_diff_z.png"))
    pylab.clf()


class IntensityAnalyser(object):
  ''' Analyse the intensities. '''

  def __init__(self, directory):
    ''' Set up the directory. '''
    from os.path import join

    # Set the directory
    self.directory = join(directory, "intensity")
    ensure_directory(self.directory)

    # Set the required fields
    self.required = [
      "intensity.raw.value",
      "intensity.raw.variance",
      "xyzcal.px",
    ]

  def __call__(self, rlist):
    ''' Analyse the reflection centroids. '''

    # Check we have the required fields
    print "Analysing reflection intensities"
    if not ensure_required(rlist, self.required):
      return

    selection = rlist['intensity.raw.variance'] < 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print 'Removing %d reflections with negative variance' % \
        selection.count(True)

    selection = rlist['intensity.raw.value'] < 0
    if selection.count(True) > 0:
      rlist.del_selected(selection)
      print 'Removing %d reflections with negative intensity' % \
        selection.count(True)

    # Select only integrated reflections
    Command.start(" Selecting only integated reflections")
    mask = rlist.get_flags(rlist.flags.integrated)
    rlist = rlist.select(mask)
    Command.end(" Selected %d integrated reflections" % len(rlist))

    # Look at distribution of I/Sigma
    print " Analysing distribution of I/Sigma"
    self.i_over_s_hist(rlist)
    print " Analysing distribution of I/Sigma vs xy"
    self.i_over_s_vs_xy(rlist)
    print " Analysing distribution of I/Sigma vs z"
    self.i_over_s_vs_z(rlist)
    print " Analysing number of background pixels used"
    self.num_background_hist(rlist)
    print " Analysing number of foreground pixels used"
    self.num_foreground_hist(rlist)

  def i_over_s_hist(self, rlist):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    pylab.title("Log I/Sigma histogram")
    pylab.hist(flex.log(I_over_S), bins=20)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "ioversigma_hist"))
    pylab.clf()

  def i_over_s_vs_xy(self, rlist):
    ''' Plot I/Sigma vs X/Y '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of I/Sigma vs X/Y")
    cax = pylab.hexbin(x, y, C=flex.log(I_over_S), gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("Log I/Sigma")
    pylab.savefig(join(self.directory, "ioversigma_vs_xy.png"))
    pylab.clf()

  def i_over_s_vs_z(self, rlist):
    ''' Plot I/Sigma vs Z. '''
    from dials.array_family import flex
    from os.path import join
    I = rlist['intensity.raw.value']
    I_sig = flex.sqrt(rlist['intensity.raw.variance'])
    I_over_S = I / I_sig
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Distribution of I/Sigma vs X/Y")
    cax = pylab.hexbin(z, flex.log(I_over_S), gridsize=100)
    pylab.xlabel("z")
    pylab.ylabel("Log I/Sigma")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "ioversigma_vs_z.png"))
    pylab.clf()

  def num_background_hist(self, rlist):
    ''' Analyse the number of background pixels. '''
    from dials.array_family import flex
    from os.path import join
    N = rlist['n_background']
    pylab.title("Num Background Pixel Histogram")
    pylab.hist(N, bins=20)
    pylab.xlabel("Number of pixels")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "n_background_hist.png"))
    pylab.clf()

  def num_foreground_hist(self, rlist):
    ''' Analyse the number of foreground pixels. '''
    from dials.array_family import flex
    from os.path import join
    N = rlist['n_foreground']
    pylab.title("Num Foreground Pixel Histogram")
    pylab.hist(N, bins=20)
    pylab.xlabel("Number of pixels")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "n_foreground_hist.png"))
    pylab.clf()


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

    # Analyse distribution of reference spots
    print " Analysing reference profile distribution vs x/y"
    self.reference_xy(rlist)
    print " Analysing reference profile distribution vs z"
    self.reference_z(rlist)

    # Look at correlations between profiles
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
    cax = pylab.hexbin(x, y, gridsize=100)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar = pylab.colorbar(cax)
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "reference_xy.png"))
    pylab.clf()

  def reference_z(self, rlist):
    ''' Analyse the distribution of reference profiles. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reference profiles binned in Z")
    cax = pylab.hist(z, bins=20)
    pylab.xlabel("z")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "reference_z.png"))
    pylab.clf()

  def reflection_corr_hist(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    pylab.title("Reflection correlations histogram")
    pylab.hist(corr, bins=20)
    pylab.xlabel("Correlation with reference profile")
    pylab.ylabel("# reflections")
    pylab.savefig(join(self.directory, "%s_corr_hist" % filename))
    pylab.clf()

  def reflection_corr_vs_xy(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations binned in X/Y")
    cax = pylab.hexbin(x, y, C=corr, gridsize=100, vmin=0.0, vmax=1.0)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("x")
    pylab.ylabel("y")
    cbar.ax.set_ylabel("Correlation with reference profile")
    pylab.savefig(join(self.directory, "%s_corr_vs_xy.png" % filename))
    pylab.clf()

  def reflection_corr_vs_z(self, rlist, filename):
    ''' Analyse the correlations. '''
    from dials.array_family import flex
    from os.path import join
    corr = rlist['profile.correlation']
    x, y, z = rlist['xyzcal.px'].parts()
    pylab.title("Reflection correlations vs Z")
    cax = pylab.hexbin(z, corr, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("z")
    pylab.ylabel("Correlation with reference profile")
    cbar.ax.set_ylabel("# reflections")
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
    cax = pylab.hexbin(flex.log(I_over_S), corr, gridsize=100)
    cbar = pylab.colorbar(cax)
    pylab.xlabel("Log I/Sigma")
    pylab.ylabel("Correlation with reference profile")
    cbar.ax.set_ylabel("# reflections")
    pylab.savefig(join(self.directory, "%s_corr_vs_ios.png" % filename))
    pylab.clf()


class Analyser(object):
  ''' Helper class to do all the analysis. '''

  def __init__(self, directory):
    ''' Setup the analysers. '''
    from os.path import join
    directory = join(directory, "analysis")
    self.analysers = [
      CentroidAnalyser(directory),
      IntensityAnalyser(directory),
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
