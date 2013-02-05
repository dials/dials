

from dials.util import MyMap, MyData

m = MyMap()
m.insert(1, MyData(1, 0))
m.insert(1, MyData(1, 0))
m.insert(1, MyData(1, 0))
m.insert(1, MyData(1, 0))
m.insert(1, MyData(1, 0))
m.insert(2, MyData(1, 0))
m.insert(2, MyData(1, 0))
m.insert(2, MyData(1, 0))

for kv in m:
    kv.second.value1 += 10
    kv.second.value2 = 100
    kv.second = MyData(-1, -1)
    
print "Access all elements"
for kv in m:
    print kv.first, kv.second.value1, kv.second.value2

for kv in m[1]:
    kv.second.value1 *= 1000

m.erase(2)

print "Access all elements"
for kv in m:
    print kv.first, kv.second.value1, kv.second.value2

#print "Access elements with index 1"
#for i, (k, v) in enumerate(m[1]):
#    m[1][i] = v.value2 + 1
#    print k, v.value1, v.value2

#print "Access all elements"
#for k, v in m:
#    print k, v.value1, v.value2


    

