
import cPickle as pickle
#v1, c1, l1 = pickle.load(open('/home/upc86896/Data/TRP_M1S3_2_/temp.p', 'rb'))
v2, c2, l2 = pickle.load(open('/home/upc86896/Data/TRP_M1S3_2_/temp.p', 'rb'))

#c3 = c1.select(l1 == 120769)
#c4 = c2.select(l2 == 120769)

#print len(c3), len(c4)

#z1, y1, x1 = zip(*c3)
#minx1, maxx1 = min(x1), max(x1)
#miny1, maxy1 = min(y1), max(y1)
#minz1, maxz1 = min(z1), max(z1)
#size1 = (maxz1 - minz1 + 1, maxy1 - miny1 + 1, maxx1 - minx1 + 1)

#z2, y2, x2 = zip(*c4)
#minx2, maxx2 = min(x2), max(x2)
#miny2, maxy2 = min(y2), max(y2)
#minz2, maxz2 = min(z2), max(z2)
#size2 = (maxz2 - minz2 + 1, maxy2 - miny2 + 1, maxx2 - minx2 + 1)

#from scitbx.array_family import flex
#a1 = flex.int(flex.grid(size1), 0)
#a2 = flex.int(flex.grid(size2), 0)

#for c in c3:
#    z = c[0] - minz1
#    y = c[1] - miny1
#    x = c[2] - minx1
#    a1[z, y, x] = 1
#
#for c in c4:
#    z = c[0] - minz2
#    y = c[1] - miny2
#    x = c[2] - minx2
#    a2[z, y, x] = 1
#    if z == 10:
#      print c
#
#a1 = a1.as_numpy_array()
#a2 = a2.as_numpy_array()

#print a1
#print a2
#print a1.shape, a2.shape

from dials.model.data import PixelList
pl = PixelList((2527, 2463), (0, 500), v2, c2)
print max(pl.labels_3d())
