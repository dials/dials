from boost.python import streambuf
from dxtbx import read_uint16
import sys
from dxtbx.format.Registry import Registry
    
format = Registry.find(sys.argv[1])
i = format(sys.argv[1])
size = i.get_detector().get_image_size()

f = open(sys.argv[1], 'rb')
hdr = f.read(512)
l = read_uint16(streambuf(f), int(size[0] * size[1]))
print sum(l)
