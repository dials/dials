

if __name__ == '__main__':
  import sys
  from dials.model.serialize import load
  from matplotlib import pylab
  rlist = load.reflections(sys.argv[1])
  I = [r.intensity for r in rlist]
  pylab.plot(I)
  pylab.show()
