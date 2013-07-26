from __future__ import division
from dials.scratch.luiso_s import add_2d, write_2d, subtrac_bkg_2d
from scitbx.array_family import flex
big_nrow = 3
big_ncol = 3

print "big (nrow, ncol) =", big_nrow, big_ncol
big_nrow = big_nrow * 2 + 1
big_ncol = big_ncol * 2 + 1
sumation = flex.double(flex.grid(big_nrow, big_ncol))
print "summation  ="
#write_2d(sumation)
descr = flex.double(flex.grid(1, 3))
'''
for ref in reflections:
    if ref.is_valid():
        
shoebox = ref.shoebox
mask = ref.shoebox_mask
background = ref.shoebox_background

data2d = shoebox[0:1, :, :]
mask2d = mask[0:1, :, :]
background2d = background[0:1, :, :]
'''
data2d = flex.double(flex.grid(3, 3))
mask2d = flex.double(flex.grid(3, 3))
background2d = flex.double(flex.grid(3, 3))

data2d[1, 1] = 1
'''
double centr_col = descriptor(0,0);
double centr_row = descriptor(0,1);
double scale = descriptor(0,2);
'''
descr[0, 0] = 1.5
descr[0, 1] = 1.666666

descr[0, 2] = 1.0
print "background2d ="
#write_2d(background2d)
peak2d = subtrac_bkg_2d(data2d, background2d)
print "peak 2d ="
write_2d(peak2d)

sumation = add_2d(descr, peak2d, sumation)
print "summation (after) ="
#from matplotlib import pyplot as plt
#print "Plotting reslt"
#img_suma = sumation.as_numpy_array()
#plt.imshow(img_suma, interpolation = "nearest")
#plt.show()
write_2d(sumation)
print "_____________________________________________________________________________________________"



print "hi 02"
