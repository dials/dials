from __future__ import division
from libtbx.phil import parse
import sys
import os
arg_defaults = '''
resolutionizer {
  rmerge = 0.0
    .type = float
  completeness = 0.0
    .type = float
  isigma = 1.0
    .type = float
  misigma = 2.0
    .type = float
  nbins = 100
    .type = int
}
'''
var1_phil = parse(arg_defaults)
arg_lst = arg_lst = sys.argv[1:]

for arg in arg_lst[1:]:
    if os.path.exists(arg):
        print "here"
        var1_phil = var1_phil.fetch(
            source = parse(open(arg).read()))
    else:
        arg_lst.append(arg)
        print "here else"

print arg_lst

#for phil_arg in phil_args:
#    interp = var1_phil.command_line_argument_interpreter(
#        home_scope = 'resolutionizer')
#    more_phil = interp.process(phil_arg)
#    var1_phil = var1_phil.fetch(source = more_phil)
#
#self._params = var1_phil.extract().resolutionizer
#
## look for useful information in the input file
#mtz_file = arg_lst[0]
#mtz_obj = mtz.object(mtz_file)
