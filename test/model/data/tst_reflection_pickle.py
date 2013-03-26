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

def tst_reflection():
    '''Test pickling the beam object.'''
    obj1 = Reflection()
    obj2 = pickle_then_unpickle(obj1)
    print "OK"
    
def tst_reflection_list():
    '''Test pickling the beam object.'''
    obj1 = ReflectionList()
    obj1.append(Reflection())
    obj1.append(Reflection())
    obj2 = pickle_then_unpickle(obj1)
    assert(len(obj2) == 2)
    print "OK"

def run():
    '''Run all the tests.'''
    tst_reflection()
    tst_reflection_list()

if __name__ == '__main__':
    run()
