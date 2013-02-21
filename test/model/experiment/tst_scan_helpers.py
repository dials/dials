
def tst_is_angle_in_range():
    """Test that for a range of angles and angular ranges, the
    is_angle_in_range function correctly calculates if the angle
    is in the range.

    """
    from dials.model.experiment import is_angle_in_range
    from random import random

    # Some helper lambda functions
    mod_360 = lambda x: x - 360.0 * floor(x / 360.0)
    mod_360_range = lambda x: (mod_360(x[0]), mod_360(x[1]))
    random_0_360 = lambda: int(random() * 360.0)

    # Create a number of ranges between 0 and 360 and check 360 degrees worth
    # of angles to see if they are within the range
    num_range = 100
    for n in range(num_range):
      angular_range = (random_0_360(), random_0_360())

      # If A < B or A > B
      if angular_range[0] < angular_range[1]:

          # Check that the following are true
          #   angle in range 0 -> A = False
          #   angle in range A -> B = True
          #   angle in range B -> 360 = False
          for angle in range(0, angular_range[0]):
              assert(is_angle_in_range(angular_range, angle, True) == False)
          for angle in range(angular_range[0], angular_range[1]+1):
              assert(is_angle_in_range(angular_range, angle, True) == True)
          for angle in range(angular_range[1]+1, 360):
              assert(is_angle_in_range(angular_range, angle, True) == False)
      else:

          # Check that the following are true
          #   angle in range 0 -> B = True
          #   angle in range B -> A = False
          #   angle in range A -> 360 = True
          for angle in range(0, angular_range[1]+1):
              assert(is_angle_in_range(angular_range, angle, True) == True)
          for angle in range(angular_range[1]+1, angular_range[0]):
              assert(is_angle_in_range(angular_range, angle, True) == False)
          for angle in range(angular_range[0], 360):
              assert(is_angle_in_range(angular_range, angle, True) == True)

    # Create a range over 360 and make sure all angles are valid
    angular_range = (-10, 370)
    for angle in range(0, 360):
        assert(is_angle_in_range(angular_range, angle, True) == True)

    # Test passed
    print "OK"

def tst_is_scan_angle_valid():
    from dials.model.experiment import Scan, is_scan_angle_valid

    # Create a scan object with no range (all angles valid)
    scan = Scan((0, 0), 0, 0.1)
    for angle in range(0, 360):
        assert(is_scan_angle_valid(scan, angle, True) == True)

    # Test passed
    print "OK"

def run():
    """Run tests for the scan_helpers.h module."""
    tst_is_angle_in_range()
    tst_is_scan_angle_valid()

if __name__ == '__main__':
    run()
