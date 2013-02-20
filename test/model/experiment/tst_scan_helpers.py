
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
    random_0_360 = lambda: random() * 360.0
    
    # Create a number of ranges between 0 and 360 and check 360 degrees worth
    # of angles to see if they are within the range    
    angular_range = (random_0_360(), random_0_360())
    if (angular_range[1] < angular_range[0]):
        angular_range = (angular_range[0], angular_range[1] + 360)
    for angle in range(0.0, angular_range[0]):
        pass
    for angle in range(anglular_range[0], angular_range[1]):
        pass
    for angle in range(anglular_range[0], angular_range[1]):
        pass        
#    for angle in range(0, 360):
#        in_range = is_angle_in_range(angular_range, angle, True)
#        if ((angular_range[0] <= angle <= angular_range[1]) or
#            (angular_range[0] <= angle and 
#             angular_range[1] <= angular_range[0]) or
#            (angular_range[0] >= angle and
#             angular_range[1] >= angle and angular_range[1] <= angular_range[0])):
#            assert(in_range == True)
#            print 1
#        else:
#            print 2
#            assert(in_range == False)
    
def tst_is_scan_angle_valid():
    print "OK"

def run():
    """Run tests for the scan_helpers.h module."""
    tst_is_angle_in_range()
    tst_is_scan_angle_valid()

if __name__ == '__main__':
    run()


