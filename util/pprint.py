

from __future__ import absolute_import, division, print_function

def profile2d(p, vmin=None, vmax=None):
  from dials.array_family import flex
  import string
  if vmin is None:
    vmin = flex.min(p)
  if vmax is None:
    vmax = flex.max(p)
  assert(vmax >= vmin)
  dv = vmax - vmin
  if dv == 0:
    c = 0
    m = 0
  else:
    m = 35.0 / dv
    c = -m * vmin
  lookup = string.digits + string.ascii_uppercase
  ny, nx = p.all()
  text = ''
  for j in range(ny):
    for i in range(nx):
      v = int(m * p[j,i] + c)
      if v < 0:
        v = 0
      elif v > 35:
        v = 35
      t = lookup[v]
      text += t + ' '
    text += '\n'
  return text

def profile3d(p, vmin=None, vmax=None):
  ''' Print a 3D profile. '''
  from dials.array_family import flex
  if vmin is None:
    vmin = flex.min(p)
  if vmax is None:
    vmax = flex.max(p)
  nz, ny, nx = p.all()
  text = []
  for k in range(nz):
    p2 = p[k:k+1,:,:]
    p2.reshape(flex.grid(ny, nx))
    text.append(profile2d(p2, vmin=vmin, vmax=vmax))
  return '\n'.join(text)

if __name__ == '__main__':

  from dials.array_family import flex
  a1 = flex.double(
    [[0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0]])
  a2 = flex.double(
    [[0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [0, 1, 2, 1, 0],
     [0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0]])
  a3 = flex.double(
    [[0, 0, 1, 0, 0],
     [0, 1, 2, 1, 0],
     [1, 2, 5, 2, 1],
     [0, 1, 2, 1, 0],
     [0, 0, 1, 0, 0]])
  a4 = flex.double(
    [[0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [0, 1, 2, 1, 0],
     [0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0]])
  a5 = flex.double(
    [[0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0]])
  a1.reshape(flex.grid(1, 5, 5))
  a2.reshape(flex.grid(1, 5, 5))
  a3.reshape(flex.grid(1, 5, 5))
  a4.reshape(flex.grid(1, 5, 5))
  a5.reshape(flex.grid(1, 5, 5))
  a = flex.double(flex.grid(5, 5, 5))
  a[0:1,:,:] = a1
  a[1:2,:,:] = a2
  a[2:3,:,:] = a3
  a[3:4,:,:] = a4
  a[4:5,:,:] = a5
  a = a * 1000

  print(profile3d(a))
