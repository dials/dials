
#refl[index[10000]].corrected_intensity
#Out[15]: 1910.0948427154015

#In [16]: refl[index[10000]].miller_index
#Out[16]: (17, 23, 23)

#In [17]: refl[index[10000]].image_coord_px
#Out[17]: (788.1612634483189, 1969.9381446581392)

#In [18]: refl[index[10000]].frame_number
#Out[18]: 352.1425769902628

#In [19]: 2.688E+03
#Out[19]: 2688.0

from dials.model.serialize import load
sweep = load.sweep('/home/upc86896/Data/TRP_M1S3_2_/sweep.json')
crystal = load.crystal('/home/upc86896/Data/TRP_M1S3_2_/crystal.json')
reference= load.reference('/home/upc86896/Data/TRP_M1S3_2_/reference.pickle')
#from dials.algorithms.integration import ReflectionExtractor

#extract = ReflectionExtractor(3)
#refl = extract(sweep, crystal)
#index = 104786
#print refl[index].miller_index

import pickle
#pickle.dump(refl[index], open('test_reflection', 'w'))
refl = pickle.load(open('test_reflection', 'r'))

from dials.model.data import ReflectionList
ref_list = ReflectionList(1)
ref_list[0] = refl

from dials.algorithms.background import XdsSubtractor
background = XdsSubtractor(min_data = 10, n_sigma = 3)
background(sweep, crystal, ref_list)

r = ref_list[0]
print r.bounding_box

shoebox = r.shoebox
background = r.shoebox_background
mask = r.shoebox_mask

diff = shoebox - background
from matplotlib import pylab
from scitbx.array_family import flex
import numpy
#for k in range(shoebox.all()[0]):
#    print diff.as_numpy_array()[k].astype(numpy.int32)

#    pylab.imshow(diff.as_numpy_array()[k], vmin = flex.min(diff), vmax = flex.max(diff))
#    pylab.show()

from dials.algorithms.integration import Summation3d, ProfileFittingReciprocalSpace
#integrate = Summation3d()
#integrate(sweep, crystal, ref_list)

from dials.algorithms.reflection_basis import transform as rbt
from dials.algorithms.integration import ProfileFittingReciprocalSpaceAlgorithm

spec = rbt.TransformSpec(sweep, crystal, 3, 4)
rbt.forward_batch(spec, ref_list)

prof = ref_list[0].transformed_shoebox
bprof = ref_list[0].transformed_shoebox_background

for k in range(9):
#    pylab.imshow(bprof.as_numpy_array()[k], vmin=flex.min(bprof), vmax=flex.max(bprof))
#    pylab.show()
    print bprof.as_numpy_array()[k].astype(numpy.int32)


print flex.sum(prof), flex.sum(shoebox)
print flex.sum(bprof), flex.sum(background)

r = ref_list[0]
coord = r.image_coord_px + (r.frame_number,)
rprof = reference.profile(coord)
print flex.sum(rprof), flex.sum(prof), flex.sum(prof - bprof)

#for k in range(9):
#    #print prof.as_numpy_array()[k].astype(numpy.int32)
#    pylab.imshow(prof.as_numpy_array()[k], vmin=flex.min(prof), vmax=flex.max(prof))
#    pylab.show()
#print ""

#for k in range(9):
##    print (3000 * rprof.as_numpy_array()[k]).astype(numpy.int32)
#    pylab.imshow(rprof.as_numpy_array()[k], vmin=flex.min(rprof), vmax=flex.max(rprof))
#    pylab.show()

# Calculate the correlation
prof_a = prof
prof_b = rprof
n = len(prof_a)
mv_a = flex.mean_and_variance(prof_a.as_1d())
mv_b = flex.mean_and_variance(prof_b.as_1d())
ma, sa = mv_a.mean(), mv_a.unweighted_sample_standard_deviation()
mb, sb = mv_b.mean(), mv_b.unweighted_sample_standard_deviation()
R = (1.0/(n-1.0)) * flex.sum((prof_a-ma) * (prof_b-mb) / (sa*sb))
print "Correlation: ", R

from dials.algorithms.reflection_basis.transform import ideal_profile
iprof = ideal_profile(4, 5)

prof_a = prof
prof_b = iprof
n = len(prof_a)
mv_a = flex.mean_and_variance(prof_a.as_1d())
mv_b = flex.mean_and_variance(prof_b.as_1d())
ma, sa = mv_a.mean(), mv_a.unweighted_sample_standard_deviation()
mb, sb = mv_b.mean(), mv_b.unweighted_sample_standard_deviation()
R = (1.0/(n-1.0)) * flex.sum((prof_a-ma) * (prof_b-mb) / (sa*sb))
print "Correlation: ", R

for k in range(9):
    pylab.subplot(3, 3, k+1)
    pylab.imshow(iprof.as_numpy_array()[k], vmin=flex.min(iprof), vmax=flex.max(iprof))
pylab.show()

print flex.sum(iprof)

rprof = iprof

from dials.algorithms.integration.profile import ProfileFitting
fit = ProfileFitting(rprof, prof, bprof, 10)
assert(fit.niter() < 10)
ref_list[0].intensity = fit.intensity()
ref_list[0].variance = fit.variance()
#ref_list[0].intensity = flex.sum(prof) - flex.sum(bprof)

#integrate = ProfileFittingReciprocalSpaceAlgorithm(reference)
#integrate(ref_list)


from dials.algorithms.integration.lp_correction import correct_intensity
correct_intensity(sweep, crystal, ref_list)
print "Intensity: ", ref_list[0].intensity
print "Corrected Intensity: ", ref_list[0].corrected_intensity
print "XDS Intensity: ", 2688.0
