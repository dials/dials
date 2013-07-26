from __future__ import division
from dials.model.data import Reflection, ReflectionList

def pickle_then_unpickle(obj):
    '''Pickle to a temp file then un-pickle.'''
    import pickle
    import tempfile

    # Create a tmp file
    temp = tempfile.TemporaryFile()

    # Pickle the object
    pickle.dump(obj, temp)

    # Read the object
    temp.flush()
    temp.seek(0)
    return pickle.load(temp)

def generate_profile(flex_type):
    '''Generate a profile array'''
    from scitbx.array_family import flex
    from random import randint, uniform, choice
    shape = randint(0, 10), randint(1, 10), randint(1, 10)
    p = flex_type(flex.grid(shape))
    for i in range(len(p)):
        p[i] = randint(0, 10000)
    return p

def generate_reflection():
    '''Generate a reflection with random attributes'''
    from random import randint, uniform, choice
    from scitbx.array_family import flex
    from math import pi
    r = Reflection()
    r.miller_index = randint(-10, 10), randint(-10, 10), randint(-10, 10)
    r.entering = choice([True, False])
    r.status = randint(0, 10)
    r.rotation_angle = uniform(-pi, pi)
    r.beam_vector = uniform(-10, 10), uniform(-10, 10), uniform(-10, 10)
    r.image_coord_mm = uniform(0, 200), uniform(0, 200)
    r.image_coord_px = uniform(0, 2000), uniform(0, 2000)
    r.frame_number = uniform(0, 1000)
    r.panel_number = randint(0, 100)
    r.centroid_position = uniform(0, 200), uniform(0, 200), uniform(0, 200)
    r.centroid_variance = uniform(0, 10), uniform(0, 10), uniform(0, 10)
    r.centroid_sq_width = uniform(0, 10), uniform(0, 10), uniform(0, 10)
    r.intensity = uniform(0, 10000)
    r.intensity_variance = uniform(0, 10000)
    r.shoebox = generate_profile(flex.double)
    r.shoebox_mask = generate_profile(flex.int)
    r.shoebox_background = generate_profile(flex.double)
    r.transformed_shoebox = generate_profile(flex.double)
    return r

def check_profile(p1, p2):
    assert(p1.all() == p2.all())
    if len(p1) > 0:
        for v1, v2 in zip(p1, p2):
            print v1, v2
            assert(v1 == v2)

def tst_reflection():
    '''Test pickling the beam object.'''
    obj1 = Reflection()
    obj2 = pickle_then_unpickle(obj1)
    print "OK"

def tst_reflection_list():
    '''Test pickling the beam object.'''
    obj1 = ReflectionList()
    for i in range(1000):
        obj1.append(generate_reflection())
    obj2 = pickle_then_unpickle(obj1)
    assert(len(obj2) == 1000)
    for r1, r2 in zip(obj1, obj2):
        assert(r1.miller_index == r2.miller_index)
        assert(r1.entering == r2.entering)
        assert(r1.status == r2.status)
        assert(r1.rotation_angle == r2.rotation_angle)
        assert(r1.beam_vector == r2.beam_vector)
        assert(r1.image_coord_mm == r2.image_coord_mm)
        assert(r1.image_coord_px == r2.image_coord_px)
        assert(r1.frame_number == r2.frame_number)
        assert(r1.panel_number == r2.panel_number)
        assert(r1.centroid_position == r2.centroid_position)
        assert(r1.centroid_variance == r2.centroid_variance)
        assert(r1.centroid_sq_width == r2.centroid_sq_width)
        assert(r1.intensity == r2.intensity)
        assert(r1.intensity_variance == r2.intensity_variance)
        check_profile(r1.shoebox, r2.shoebox)
        check_profile(r1.shoebox_mask, r2.shoebox_mask)
        check_profile(r1.shoebox_background, r2.shoebox_background)
        check_profile(r1.transformed_shoebox, r2.transformed_shoebox)

    print "OK"

def run():
    '''Run all the tests.'''
    tst_reflection()
    tst_reflection_list()

if __name__ == '__main__':
    run()
