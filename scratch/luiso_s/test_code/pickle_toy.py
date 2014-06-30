import cPickle as pickle
from dials.array_family import flex
reflections = pickle.load(open('refl_01.pickle', 'rb'))
# Get as a reflection table
table = reflections.to_table()

table.as_pickle("integrated_reflections_table.pickle")
