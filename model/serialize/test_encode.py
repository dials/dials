
#import json
#import textwrap

#scan = {
#    "image_range" : (1, 2),
#    "oscillation" : (1, 2),
#    "exposure_time" : 1.0,
#    "epochs" : [
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#    ]
#  }

#class Scan:
#    image_range = (1, 2)
#    oscillation = (1, 2)
#    exposure_time = 1.0
#    epochs = [
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#    ]

#def scan_to_dict(obj):
#    from collections import OrderedDict
#    scan = OrderedDict()
#    scan['image_range'] = obj.image_range
#    scan['oscillation'] = obj.oscillation
#    scan['exposure_time'] = obj.exposure_time
#    scan['epochs'] = obj.epochs
#    return scan
#

#class Encoder(json.JSONEncoder):
#    def default(self, obj):
#
#        if isinstance(obj, Scan):
#            d = scan_to_dict(obj)
#            return d
#
#    def iterencode(self, obj, **kwargs):
#        return json.JSONEncoder.iterencode(self, obj, **kwargs)

#scan = Scan()

#print Encoder().encode(scan)
from dials.model.experiment import Scan, Beam, Goniometer, Detector, Panel, Crystal

filenames = 'filenames###.cbf'
scan = Scan((1, 100), (0, 1), 0)
beam = Beam((0, 0, 1), 1.0)
gonio = Goniometer()
detector = Detector(Panel("UNKNOWN", (1, 0, 0), (0, 1, 0), (0, 0, 1), (0.1, 0.1), (100, 100), (0, 100)))
crystal = Crystal((1, 0, 0), (0, 1, 0), (0, 0, 1), 10)

from dials.model.serialize import dumps, loads

print dumps(filenames=filenames, scan=scan, beam=beam, goniometer=gonio, crystal=crystal, detector=detector)



#print loads(s)

#for line in textwrap.wrap(json.dumps(scan, sort_keys=False, cls=Encoder), 80):
#    print line
