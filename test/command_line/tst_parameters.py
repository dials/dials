from __future__ import division

def run_test():
  from subprocess import call

  # Call dials.parameters
  call(["dials.parameters > parameters.txt"], shell=True)

  # Read the parameters
  with open("parameters.txt", "r") as infile:
    output = infile.read()

  # Check we have some output
  assert(len(output) != 0)
  print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run_test()
