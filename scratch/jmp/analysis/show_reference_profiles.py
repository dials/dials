from __future__ import division

def print_profile(profile):
  import numpy
  numpy.set_printoptions(precision=3, suppress=True)

  slices = []
  for k in range(9):
    text = str(profile[k])
    text = text.translate(None, "[]")
    lines = text.split("\n")
    lines = [l.strip() for l in lines]
    text = "\n".join(lines)
    slices.append(text)

#    text = ""
#    for k in range(9):
#        text += slices[k] + "\n"
#
  lines = []
  for k in range(9):
    lines.append(slices[k].split("\n"))

  space = " " * 10
  row1 = ""
  row2 = ""
  row3 = ""
  for k in range(9):
    row1 += lines[0][k] + space + lines[1][k] + space + lines[2][k] + "\n"
    row2 += lines[3][k] + space + lines[4][k] + space + lines[5][k] + "\n"
    row3 += lines[6][k] + space + lines[7][k] + space + lines[8][k] + "\n"

  text = row1 + "\n" + row2 + "\n" + row3 + "\n"

  return text

from matplotlib import pylab
from scitbx.array_family import flex
from dials.model.serialize import load
import numpy
import sys

numpy.set_printoptions(precision=3, suppress=True)

reference = load.reference(sys.argv[1])
outfile = open("profiles.txt", "w")

for nref in range(9):

  profile = reference.profile(nref).as_numpy_array()
  outfile.write("Profile {0}\n".format(nref))
  outfile.write(print_profile(profile))

  xcoord = int(reference.coord(nref)[0])
  ycoord = int(reference.coord(nref)[1])

  fig = pylab.figure(1, figsize=(6, 6))
  fig.text(.5, .95, "Position: {0}".format((xcoord, ycoord)), horizontalalignment='center')
  max_profile = flex.max(reference.profile(nref))
  for k in range(9):
    ax = pylab.subplot(3, 3, k+1)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    pylab.imshow(profile[k], interpolation='none', vmin=0, vmax=max_profile)
  #pylab.show()
  fig.savefig("reference_{0}.png".format(nref))
  fig.clear()
