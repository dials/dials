def read_template(template_name='DLS6MSN100.cbft'):
  return open(template_name, 'r').read()

def get_header_and_image(in_image):
  import binascii

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(in_image, 'rb').read()
  data_offset = data.find(start_tag)
  cbf_header = data[:data_offset]

  data_offset = data.find('_array_data.data')

  image_data = data[data_offset:]

  return cbf_header, image_data

def parse_cbf_header_tokens(cbf_header):
  header = { }
  for record in cbf_header.split('\n'):
    if not record.startswith('#'):
      continue
    record = record.strip()
    tokens = record.split()
    if len(tokens) == 2:
      header['Timestamp'] = tokens[1]
      continue
    keyword = tokens[1]
    value = ' '.join(tokens[2:])
    header[keyword] = value
  return header

def string_keep(str, keep):
  result = ''
  for s in str:
    if s in keep:
      result += s
  return result

def make_image_from_image(in_image, out_image,
                          template_name='DLS6MSN100.cbft'):
  header, image = get_header_and_image(in_image)
  header = parse_cbf_header_tokens(header)

  # static things...

  header['Compression_type'] = 'byte_offset'
  header['X_dimension'] = '2463'
  header['Y_dimension'] = '2527'

  template = read_template(template_name)

  template_map = { }

  for record in template.split('\n'):
    if not record.startswith('@'):
      continue
    keyword, value = record.split()[1:3]
    assert(keyword in header)
    template_map[value] = header[keyword].split()[0]

  for token in sorted(template_map):
    template = template.replace(token, template_map[token])

  template = '###CBF: VERSION 1.5, CBFlib v0.7.8 - PILATUS detectors\r\n' + \
             '\r\ndata_block\r\n' + template.split('--- End of preamble')[1]

  beam = map(float, string_keep(header['Beam_xy'], '0123456789. ').split())

  open(out_image, 'wb').write(template % {
      'distance':1000 * float(header['Detector_distance'].split()[0]),
      'beamline':'I04',
      'xtal_id':'xtal',
      'detector_id':'P6M_I04',
      'detector_name':'Pilatus 6M',
      'beam_x':beam[0] * 0.172,
      'beam_y':beam[1] * 0.172,
      'pixel_x':0.172,
      'pixel_y':0.172
      } + image)

  return

if __name__ == '__main__':
  import sys
  make_image_from_image(sys.argv[1], sys.argv[2])
