

fullpath = '/home/upc86896/Data/TRP_M1S3_2_/reference.pickle'
small_template = '/home/upc86896/Data/TRP_M1S3_2_/reference%d.pickle'
smallpaths = [small_template % i for i in range(5)]

from dials.model.serialize import load

full_reference = load.reference(fullpath)
small_reference = [load.reference(path) for path in smallpaths]

profile1 = full_reference.profile(0)
profile2 = small_reference[0].profile(0)

from scitbx.array_family import flex
diff = flex.abs(profile1 - profile2)
print list(diff)
assert(diff.all_lt(1e-5))
