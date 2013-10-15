

from iotbx.xds import integrate_hkl
import sys

reader = integrate_hkl.reader()
reader.read_file(sys.argv[1])

hkl = []
iobs = []
sigma = []
xyzcal = []
rlp = []
peak = []
corr = []
maxc = []
xyzobs = []
alfbet0 = []
alfbet1 = []
psi = []
for i in range(len(reader.xyzcal)):
    x, y, z = reader.xyzcal[i]
    if z >= 0 and z < 9:
        hkl.append(reader.hkl[i])
        iobs.append(reader.iobs[i])
        sigma.append(reader.sigma[i])
        xyzcal.append(reader.xyzcal[i])
        rlp.append(reader.rlp[i])
        peak.append(reader.peak[i])
        corr.append(reader.corr[i])
        maxc.append(reader.maxc[i])
        xyzobs.append(reader.xyzobs[i])
        alfbet0.append(reader.alfbet0[i])
        alfbet1.append(reader.alfbet1[i])
        psi.append(reader.psi[i])

print "Num: ", len(hkl)

f = open(sys.argv[1], 'r')
lines = f.readlines()
del f
header = []
for l in lines:
    if l.strip().startswith('!END_OF_HEADER'):
        header.append('!END_OF_HEADER')
        break
    else:
        header.append(l)
header_text = ''.join(header)
all_text = header_text
for i in range(len(hkl)):
    line = list(hkl[i]) + [iobs[i]] + [sigma[i]] + list(xyzcal[i]) + [rlp[i]] + \
        [peak[i]] + [corr[i]] + [maxc[i]] + list(xyzobs[i]) + \
        list(alfbet0[i]) + list(alfbet1[i]) + [psi[i]]
    text = ' '.join([str(l) for l in line])
    all_text += text + '\n'
all_text += '!END_OF_DATA'

f = open(sys.argv[2], 'w')
f.write(all_text)
