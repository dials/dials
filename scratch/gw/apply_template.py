def find_number(words):
    for token in words.split():
        try:
            return float(token)
        except ValueError, e:
            pass

def keep_numbers(words):
    keep = ' 0123456789.-'
    result = ''
    for letter in words:
        if letter in keep:
            result += letter
    return result

template = 'i19-1.cbftt'

keys = []
name_hash = { }

for record in open(template):
    if record.startswith('@'):
        tokens = record.split()
        keys.append(tokens[1])
        name_hash[tokens[2]] = tokens[1]

import sys

cbf_file = sys.argv[1]

prior = {'X_dimension':1475,
         'Y_dimension':1679,
         'Compression_type':'CBF_BYTE_OFFSET'}

found = { }

for key in prior:
    found[key] = prior[key]

beam_x = None
beam_y = None

for record in open(cbf_file):
    if '_array_data.data' in record:
        break
    if record.startswith('#'):
        found[record.split()[1].replace(':', '')] = find_number(record)
        if len(record.replace('#', '').strip().split()) == 1:
            found['Timestamp'] = record.replace('#', '').strip().split()[0]
        if 'Beam_xy' in record:
            beam_x, beam_y = map(float, keep_numbers(record).split())

template_str = open(template).read()

import os

for key in found:
    if key in os.environ:
        found[key] = os.environ[key]

for key in name_hash:
    template_str = template_str.replace(key, str(found[name_hash[key]]))

# scale beam, hard code default

if beam_x > 0:
    beam_x *= 0.172
else:
    beam_x = 730.00 * 0.172

if beam_y > 0:
    beam_y *= 0.172
else:
    beam_y = 865.00 * 0.172

values = {'distance':1000*found['Detector_distance'],
          'pixel_x':0.172,
          'pixel_y':0.172,
          'beamline':'I19-1',
          'detector_id':'P2M',
          'detector_name':'P2M S/N 24-0107',
          'xtal_id':'XTAL',
          'beam_x':beam_x,
          'beam_y':beam_y}

template_text = (template_str % values).split('End of preamble')[-1].strip()

raw_data = open(cbf_file, 'rb').read()
header = raw_data.split('_array_data.data')[0]
data = raw_data.split('_array_data.data')[-1]

new_data = header + template_text + '\r\n_array_data.data' + data

assert not os.path.exists(sys.argv[2])
open(sys.argv[2], 'wb').write(new_data)
