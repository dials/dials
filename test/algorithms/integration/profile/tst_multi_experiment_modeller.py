
from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    from dials.algorithms.profile_model.modeller import ProfileModellerIface
    from dials.algorithms.profile_model.modeller import MultiExpProfileModeller
    from dials.array_family import flex

    class Modeller(ProfileModellerIface):

      def __init__(self, index, expected):
        self.index = index
        self.accumulated = False
        self.finalized = False
        self.expected = expected
        super(Modeller, self).__init__()

      def model(self, reflections):
        assert(reflections['id'].all_eq(self.index))
        assert(len(reflections) == self.expected)

      def accumulate(self, other):
        self.accumulated = True
        assert(self.index == other.index)

      def finalize(self):
        assert(self.accumulated)
        self.finalized = True

    # The expected number of reflections
    expected = [100, 200, 300, 400, 500]

    # Create some reflections
    reflections = flex.reflection_table()
    reflections["id"] = flex.int()
    for idx in range(len(expected)):
      for n in range(expected[idx]):
        reflections.append({
          "id" : idx
        })

    # Create two modellers
    modeller1 = MultiExpProfileModeller()
    modeller2 = MultiExpProfileModeller()
    for idx in range(len(expected)):
      modeller1.add(Modeller(idx, expected[idx]))
      modeller2.add(Modeller(idx, expected[idx]))

    # Model the reflections
    modeller1.model(reflections)
    modeller2.model(reflections)

    # Accumulate
    modeller1.accumulate(modeller2)

    # Finalize
    modeller1.finalize()

    # Check finalized
    assert(modeller1.finalized)

    # Test passed
    print 'OK'


if __name__ == '__main__':

  test = Test()
  test.run()
