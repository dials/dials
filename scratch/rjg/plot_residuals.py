from __future__ import division
from cctbx.array_family import flex
import math

def run(args):
  centroids_filename = args[0]
  hkl = flex.miller_index()
  frame_obs = flex.double()
  x_obs = flex.double()
  y_obs = flex.double()
  phi_obs = flex.double()
  x_calc = flex.double()
  y_calc = flex.double()
  phi_calc = flex.double()
  with open(centroids_filename, 'rb') as f:
    for i, line in enumerate(f.readlines()):
      tokens = line.split()
      if i == 0:
        print tokens
        assert tokens == [
          'H','K', 'L', 'Frame_obs', 'X_obs', 'Y_obs', 'Phi_obs', 'X_calc',
          'Y_calc', 'Phi_calc']
      else:
        hkl.append([int(t) for t in tokens[:3]])
        frame_obs.append(float(tokens[3]))
        x_obs.append(float(tokens[4]))
        y_obs.append(float(tokens[5]))
        phi_obs.append(float(tokens[6]))
        x_calc.append(float(tokens[7]))
        y_calc.append(float(tokens[8]))
        phi_calc.append(float(tokens[9]))

  phi_obs_deg = (180/math.pi) * phi_obs

  x_residuals = x_calc - x_obs
  y_residuals = y_calc - y_obs
  phi_residuals = phi_calc - phi_obs

  mean_residuals_x = []
  mean_residuals_y = []
  mean_residuals_phi = []
  phi = []

  for i_phi in range(int(math.floor(flex.min(phi_obs_deg))),
                 int(math.ceil(flex.max(phi_obs_deg)))):
    sel = (phi_obs_deg >= i_phi) & (phi_obs_deg < (i_phi+1))
    if sel.count(True) == 0:
      continue
    mean_residuals_x.append(flex.mean(x_residuals.select(sel)))
    mean_residuals_y.append(flex.mean(y_residuals.select(sel)))
    mean_residuals_phi.append(flex.mean(phi_residuals.select(sel)))
    phi.append(i_phi)

  from matplotlib import pyplot
  fig = pyplot.figure()
  ax = fig.add_subplot(311)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_x)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_x')
  ax = fig.add_subplot(312)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_y)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_y')
  ax = fig.add_subplot(313)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_phi)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_phi')
  pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
