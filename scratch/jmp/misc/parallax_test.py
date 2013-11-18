

if __name__ == '__main__':
#    from scitbx import matrix
#    from scitbx.array_family import flex
#    from iotbx.xds import correction
#    from cctbx import factor_ev_angstrom, factor_kev_angstrom

  factor_mev_angstrom = (6.6260755 * 2.99792458 / 1.60217733) / 1000.0

  width, height = (2463, 2527)
  wavelength = 0.979500

#    print "Reading X correction"
#    handle = correction.reader()
#    handle.read_file('/home/upc86896/Data/X4_lots_M1S4_1_/GX-CORRECTIONS.cbf')
#    xcorr = handle.get_correction((height, width))
#
#    print "Reading Y correction"
#    handle = correction.reader()
#    handle.read_file('/home/upc86896/Data/X4_lots_M1S4_1_/GY-CORRECTIONS.cbf')
#    ycorr = handle.get_correction((height, width))


  table = [[1.00000E-03,  1.570E+03,  1.567E+03],
            [1.50000E-03,  5.355E+02,  5.331E+02],
            [1.83890E-03,  3.092E+02,  3.070E+02],
            [1.83890E-03,  3.192E+03,  3.059E+03],
            [2.00000E-03,  2.777E+03,  2.669E+03],
            [3.00000E-03,  9.784E+02,  9.516E+02],
            [4.00000E-03,  4.529E+02,  4.427E+02],
            [5.00000E-03,  2.450E+02,  2.400E+02],
            [6.00000E-03,  1.470E+02,  1.439E+02],
            [8.00000E-03,  6.468E+01,  6.313E+01],
            [1.00000E-02,  3.389E+01,  3.289E+01],
            [1.50000E-02,  1.034E+01,  9.794E+00],
            [2.00000E-02,  4.464E+00,  4.076E+00],
            [3.00000E-02,  1.436E+00,  1.164E+00],
            [4.00000E-02,  7.012E-01,  4.782E-01],
            [5.00000E-02,  4.385E-01,  2.430E-01],
            [6.00000E-02,  3.207E-01,  1.434E-01],
            [8.00000E-02,  2.228E-01,  6.896E-02],
            [1.00000E-01,  1.835E-01,  4.513E-02],
            [1.50000E-01,  1.448E-01,  3.086E-02],
            [2.00000E-01,  1.275E-01,  2.905E-02],
            [3.00000E-01,  1.082E-01,  2.932E-02],
            [4.00000E-01,  9.614E-02,  2.968E-02],
            [5.00000E-01,  8.748E-02,  2.971E-02],
            [6.00000E-01,  8.077E-02,  2.951E-02],
            [8.00000E-01,  7.082E-02,  2.875E-02],
            [1.00000E+00,  6.361E-02,  2.778E-02],
            [1.25000E+00,  5.688E-02,  2.652E-02],
            [1.50000E+00,  5.183E-02,  2.535E-02],
            [2.00000E+00,  4.480E-02,  2.345E-02],
            [3.00000E+00,  3.678E-02,  2.101E-02],
            [4.00000E+00,  3.240E-02,  1.963E-02],
            [5.00000E+00,  2.967E-02,  1.878E-02],
            [6.00000E+00,  2.788E-02,  1.827E-02],
            [8.00000E+00,  2.574E-02,  1.773E-02],
            [1.00000E+01,  2.462E-02,  1.753E-02],
            [1.50000E+01,  2.352E-02,  1.746E-02],
            [2.00000E+01,  2.338E-02,  1.757E-02]]


  from math import log, exp
  import numpy
  from scipy.interpolate import interp1d
  energy, mac, meac = zip(*table)
  energy = numpy.array(energy, dtype=numpy.float32)
  mac = numpy.array(mac, dtype=numpy.float32)

  beam_energy = factor_mev_angstrom / wavelength


#    from matplotlib import pylab
#    pylab.plot(energy, mac)
#    pylab.show()


  print energy
  print mac

#    x0 = 0
#    x1 = 0
#    x2 = 0
#
#    y0 = 0
#    y1 = 0
#    y2 = 0

#    a11 = 2.0 / (x1 - x0)
#    a12 = 1.0 / (x1 - x0)
#    a21 = 1.0 / (x1 - x0)
#    a22 = 2.0 * (1.0 / (x1 - x0) + 1.0 / (x2 - x1))
#    a23 = 1.0 / (x2 - x1)
#    a32 = 1.0 / (x2 - x1)
#    a33 = 2.0 / (x2 - x1)
#    b1 = 3.0 * (y1 - y0) / (x1 - x0)**2
#    b2 = 3.0 * ((y1 - y0) / (x1 - x0)**2 + (y2 - y1) / (x2 - x1)**2)
#    b3 = 3.0 * (y2 - y1) / (x2 - x1)**2

  log_energy = numpy.log(energy)
  log_mac = numpy.log(mac)

  print log(beam_energy)

#    from matplotlib import pylab
#    pylab.plot(energy, log_mac)
#    pylab.show()


#    print 1/0

#    from scipy.interpolate import splrep, splev
#    tck = splrep(log_energy, log_mac)
#    print tck
#    print 1/0
#    coeff = splev(log(beam_energy), tck, der=0)
#    f2 = interp1d(numpy.log(energy), numpy.log(mac), kind='linear')
#    coeff = exp(f2(log(beam_energy)))
  coeff = numpy.exp(numpy.interp(log(beam_energy), numpy.log(energy), numpy.log(mac)))
#    coeff = 16.9612
  rho = 2.330E+00

  print "Energy: {0} Mev".format(beam_energy)
  print "Energy: {0} ev".format(beam_energy * 1000000.0)
  print "Mu/rho: {0} cm^2/g".format(coeff)
  print "Mu:     {0} cm^-1".format(coeff * rho)
  print "rho:    {0} g / cm^3".format(rho)
  from math import exp
  x_arr = []
  p_arr = []


  print "Att Len: {0} microm".format((1.0 / (coeff * rho)) * 10 * 1000)

  for xx in range(0, 100):
    xcm = xx / 1000.0
    xmm = xcm * 10
    p = exp(-coeff * rho * xcm)
    x_arr.append(xmm)
    p_arr.append(p)

#    from matplotlib import pylab
#    pylab.plot(x_arr, p_arr)
#    pylab.axhline(1.0/exp(1.0))
#    pylab.show()
