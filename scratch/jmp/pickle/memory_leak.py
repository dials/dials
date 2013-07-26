

from dials.model.data import ReflectionList, flex_flex_double
from time import sleep
from psutil import virtual_memory
from scitbx.array_family import flex
import gc

def create_and_destroy():

    rlist = ReflectionList(100000)
#    for i in range(len(rlist)):
##        rlist[i].miller_index = (1, 2, 3)
#        r= rlist[i]
#        r.shoebox = flex.double(1000)

    for r in rlist:
#        r.miller_index = (1, 2, 3)
        r.shoebox = flex.double(1000)

    print virtual_memory()
    del rlist
    del r
    print virtual_memory()

print virtual_memory()
#for i in range(10):
create_and_destroy()
print virtual_memory()
for i in range(5):
    print virtual_memory()
    sleep(0.5)
