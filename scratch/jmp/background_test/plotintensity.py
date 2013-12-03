

if __name__ == '__main__':
  import sys
  from dials.model.serialize import load
  from matplotlib import pylab
  from scitbx.array_family import flex
  from dials.algorithms.shoebox import MaskCode
  from math import sqrt

  print "Loading Reflections"
  rlist = load.reflections(sys.argv[1])
  expected_intensity = float(sys.argv[2])
  expected_background = float(sys.argv[3])

  print "Calculating Errors"
#  I = [r.intensity for r in rlist]
#  E = []
#  ios = []
  I = []
  V = []
  for r in rlist:

    background_pixels = r.shoebox.select(flex.bool(
      [bool(p & MaskCode.BackgroundUsed) for p in r.shoebox_mask]))

    #background_pixels = r.shoebox.as_1d()

    signal_pixels = r.shoebox.select(flex.bool(
      [bool(p & MaskCode.Foreground) for p in r.shoebox_mask]))

    m = len(signal_pixels)
    n = len(background_pixels)
    Ib = m*flex.mean(r.shoebox_background)#(m / n) * flex.sum(background_pixels)
    Is = flex.sum(signal_pixels) - Ib
    intensity = Is
    variance = Is +Ib + (m / n)*Ib

#    pixels = r.shoebox.select(flex.bool(
#      [bool(p & MaskCode.BackgroundUsed) for p in r.shoebox_mask]))
#    mv = flex.mean_and_variance(pixels)
#    mean = mv.mean()
#    var = mv.unweighted_sample_variance()
#    sdev = sqrt(var)
#    err = sdev / sqrt(len(pixels))
#    n_background = len(pixels)
#    n_foreground = [bool(p & MaskCode.Foreground) for p in r.shoebox_mask].count(True)
#    bg_err = sqrt(n_foreground * err**2)
#    tot_err = sqrt(n_foreground * err**2 + r.intensity_variance)
#    E.append(tot_err)
#    ios.append(r.intensity / tot_err)
#    print ("Background: %f, "
#           "Intensity: %f, "
#           "Sigma: %f, "
##           "N Background: %d, "
##           "N Foreground: %d, "
#           "Background Error: %f, "
#           "Tot Background Error: %f, "
#           "Background Within Error: %s "
#           "I over Sigma: %f") % (
#      flex.mean(r.shoebox_background), r.intensity, sqrt(r.intensity_variance),
#      #n_background, n_foreground,
#      err, bg_err,
#      (flex.mean(r.shoebox_background) - expected_background) < err,
#      r.intensity / tot_err)

#    pixels = r.shoebox.select(flex.bool(
#      [bool(p & MaskCode.Foreground) for p in r.shoebox_mask]))
#    subtracted = pixels - flex.mean(r.shoebox_background)
#    print flex.mean(pixels)
    print m, n, intensity, sqrt(variance), r.intensity, sqrt(r.intensity_variance)
    I.append(intensity)
    V.append(variance)

  S = [sqrt(v) for v in V]
  Z = [(i - expected_intensity) / s for i, s in zip(I, S)]

  mv = flex.mean_and_variance(flex.double(Z))
  meanz = mv.mean()
  varz = mv.unweighted_sample_variance()
  sdevz = sqrt(varz)

  print "Z: (mean: %f), (var: %f), (sdev: %f)" % (meanz, varz, sdevz)

  pylab.hist(Z)
  pylab.show()


#  pylab.errorbar(range(len(I)), I, yerr=E, fmt='o')
#  pylab.axhline(expected_intensity)
#  pylab.show()
#
#  pylab.scatter(range(len(ios)), ios)
#  pylab.axhline(0)
#  pylab.show()
