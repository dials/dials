


from dials.model.serialize import load
import pickle

sweep = load.sweep('/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/sweep.json')

text = pickle.dumps(sweep)

sweep2 = pickle.loads(text)

sweep2.get_detector()

