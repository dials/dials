from __future__ import absolute_import, division
from __future__ import print_function
import webbrowser
import urllib2
import os

tempdir = os.getcwd()

tutorial_data = 'https://zenodo.org/record/10271/files/th_8_2.tar.bz2'
tutorial_html = 'http://dials.github.io/documentation/tutorials/processing_in_detail_tutorial.html'

print('Opening %s in web browser' % tutorial_html)
webbrowser.open(tutorial_html)

file_name = tutorial_data.split('/')[-1]
u = urllib2.urlopen(tutorial_data)
f = open(os.path.join(tempdir, file_name), 'wb')
meta = u.info()
file_size = int(meta.getheaders("Content-Length")[0])
print("Downloading: %s Bytes: %s to %s" % (file_name, file_size,
                                           os.path.join(tempdir, file_name)))

file_size_dl = 0
block_sz = 2**16

while True:
  buffer = u.read(block_sz)
  if not buffer:
    break

  file_size_dl += len(buffer)
  f.write(buffer)
  status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
  status = status + chr(8)*(len(status)+1)
  print status,
f.close()

print 'You are ready to go: tar xvfj %s' % file_name
