import cPickle as pickle
from dials.array_family import flex
reflections = pickle.load(open('integrated.pickle', 'rb'))
# Get as a reflection table
table = reflections.to_table()

table.as_pickle("integrated_reflections_table.pickle")
