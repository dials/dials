

if __name__ == '__main__':
  import sys
  from dials.model.serialize import load
  from matplotlib import pylab
  from scitbx.array_family import flex
  rlist = load.reflections(sys.argv[1])
  I = [r.intensity for r in rlist]
  for r in rlist:
    print "Background: %f, Intensity: %f" % (
      flex.mean(r.shoebox_background), r.intensity)
  pylab.plot(I)
  pylab.show()
