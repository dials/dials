


from dials.model.data import ReflectionList
from time import time
import cPickle as pickle
from scitbx.array_family import flex
import resource

print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

print "Creating reflection list"
rlist = ReflectionList(100000)
for r in rlist:
    r.miller_index = (10, 10, 10)
    r.bounding_box = (0, 10, 0, 10, 0, 10)
    r.shoebox = flex.double(flex.grid(10, 10, 10))
    r.shoebox_mask = flex.int(flex.grid(10, 10, 10))

print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

print "Pickling"
st = time()
pickle.dump(rlist, open('temp.pickle', 'wb'), pickle.HIGHEST_PROTOCOL)
print "Pickled in {0}s".format(time() - st)

print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss






print "Loading"
st = time()
rlist = pickle.load(open("temp.pickle", "rb"))
print rlist[0].miller_index
print rlist[-1].miller_index
print rlist[200].shoebox.all()
print "Loaded in {0}s".format(time() - st)
