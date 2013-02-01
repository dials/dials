import sys
from iotbx.detectors import ImageFactory

image = ImageFactory(sys.argv[1])
image.read()

nfast = image.parameters['SIZE1']
nslow = image.parameters['SIZE2']

data = image.get_raw_data()


for s in range(nslow):
    for f in range(nfast):
        print data[s * nfast + f],
    print 
