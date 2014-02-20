

def parse_frame_data(filename):
  ''' Parse the INTEGRATE.LP file. '''

  # The row headings
  header = ['IMAGE', 'IER', 'SCALE', 'NBKG', 'NOVL',
            'NEWALD', 'NSTRONG','NREJ', 'SIGMAB', 'SIGMAR']

  # Parse the file
  rows = []
  with open(filename, 'r') as infile:
    inside_table = False
    for line in infile.readlines():
      tokens = line.split()
      if len(tokens) == 10:
        if inside_table and int(tokens[0]) == len(rows) + 1:
          row = (int(tokens[0]),
                 int(tokens[1]),
                 float(tokens[2]),
                 int(tokens[3]),
                 int(tokens[4]),
                 int(tokens[5]),
                 int(tokens[6]),
                 int(tokens[7]),
                 float(tokens[8]),
                 float(tokens[9]))
          rows.append(row)
        else:
          if tokens == header:
            inside_table = True
      else:
        inside_table=False

  # Return the rows
  return header, rows


def parse_profile_data(filename):
  ''' Parse the INTEGRATE.LP file. '''

  # The row headings
  header1 = "STANDARD DEVIATION (DEGREES) OF BEAM DIVERGENCE"
  header2 = "STANDARD DEVIATION (DEGREES) OF REFLECTING RANGE"

  # Parse the file
  sigma_b_rows = []
  sigma_m_rows = []
  with open(filename, 'r') as infile:
    inside_sigmab_table = False
    inside_sigmam_table = False
    for line in infile.readlines():
      if inside_sigmab_table:
        if "BATCH" in line or "NUMBER" in line:
          continue
        tokens = map(float, line.split())
        if len(tokens) == 10:
          sigma_b_rows.append(tokens)
        else:
          inside_sigmab_table = False
      elif inside_sigmam_table:
        if "BATCH" in line or "NUMBER" in line:
          continue
        tokens = map(float, line.split())
        if len(tokens) == 10:
          sigma_m_rows.append(tokens)
        else:
          inside_sigmam_table = False
      elif header1 in line:
        inside_sigmab_table = True
      elif header2 in line:
        inside_sigmam_table = True

  # Return the rows
  return sigma_b_rows, sigma_m_rows


if __name__ == '__main__':
  import sys

  # Parse the file
  header, rows = parse_frame_data(sys.argv[1])

  from matplotlib import pylab
  sigmab = zip(*rows)[8]
  sigmam = zip(*rows)[9]

  pylab.subplot(121)
  pylab.plot(sigmab)
  pylab.subplot(122)
  pylab.plot(sigmam)
  pylab.show()

  sigma_b_rows, sigma_m_rows = parse_profile_data(sys.argv[1])

  sigma_b_data = zip(*sigma_b_rows)[1:]
  sigma_b_data = [d for e in sigma_b_data for d in e]
  sigma_m_data = zip(*sigma_m_rows)[1:]
  sigma_m_data = [d for e in sigma_m_data for d in e]
  print "Mean Sigma B: ", sum(sigma_b_data) / len(sigma_b_data)
  print "Mean Sigma M: ", sum(sigma_m_data) / len(sigma_m_data)
