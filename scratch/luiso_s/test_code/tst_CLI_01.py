from __future__ import division
import sys

arg_lst = sys.argv[1:]

algrm = 'xds'
filenames = []
times = 5
shift = 10
for arg in arg_lst:
  if '=' in arg:
    print '= in argument'
    leng_01 = arg.find('=')
    if arg[:leng_01] == 'algr' or arg[:leng_01] == 'al':
      if arg[leng_01 + 1:] == 'lui':
        algrm = 'lui'
      else:
        algrm = 'xds'
    elif arg[:leng_01] == 'times' or arg[:leng_01] == 'tm':
      times = int(arg[leng_01 + 1:])
    elif arg[:leng_01] == 'shf' or arg[:leng_01] == 'shift':
      shift = int(arg[leng_01 + 1:])
  else:
    filenames.append(arg)
print 'filenames =', filenames
print 'algrm =', algrm
print 'times=', times
print 'shift=', shift
