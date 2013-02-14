from boost.python import streambuf
from dxtbx import read_uint16
f = open('/Users/graeme/data/demo/insulin_1_001.img', 'rb')
hdr = f.read(512)
l = read_uint16(streambuf(f), 2304 * 2304)
print sum(l)
