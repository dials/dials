from __future__ import division
import numpy
t_size = 5
f_size = t_size
r_size = t_size
c_size = t_size
data3d = numpy.arange(f_size * r_size * c_size, dtype = int).reshape(f_size, r_size, c_size)
#
#print data3d
#for f in range(f_size):
#    for r in range(r_size):
#        for c in range(c_size):
#            print 'data3d[', f, ',', r, ',', c, '] =', data3d[f, r, c]
#
data = []
cont = 0
for f in range(f_size):
    for r in range(r_size):
        for c in range(c_size):
            data.append((f, r, c, cont ** 2))
            cont += 1

print 'data =', data
#f_tot = 0.0
#r_tot = 0.0
#c_tot = 0.0
#d_tot = 0.0
#
cont = 0
for f, r, c, d in data:
    print cont
    #print 'data3d[', f, ',', r, ',', c, '] =', data[f, r, c]
    print f, c, r, d
    cont += 1
#    f_tot += d * f
#    r_tot += d * r
#    c_tot += d * c
#    d_tot += d

#for f in range(f_size):
#    for r in range(r_size):
#        for c in range(c_size):
#            data3d.append((f, r, c, 0))



#for f, r, c, d in data3d:
#    print 'data3d[', f, ',', r, ',', c, '] =', d
