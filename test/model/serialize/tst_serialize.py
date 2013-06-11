
class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_crystal()

    def tst_crystal(self):
        from dials.model.serialize.crystal import crystal_to_dict
        from dials.model.serialize.crystal import crystal_from_dict
        from dials.model.experiment import Crystal

        c1 = Crystal(
            real_space_a=(1, 0, 0),
            real_space_b=(0, 1, 0),
            real_space_c=(0, 0, 1),
            sg=10,
            mosaicity=0.1)

        d = crystal_to_dict(c1)
        c2 = crystal_from_dict(d)
        assert(d['real_space_a'] == (1, 0, 0))
        assert(d['real_space_b'] == (0, 1, 0))
        assert(d['real_space_c'] == (0, 0, 1))
        assert(d['space_group'] == 10)
        assert(d['mosaicity'] == 0.1)
        assert(c1 == c2)
        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()
