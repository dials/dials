

class PartialShoeboxMerger(object):

    def __init__(self, tar, predicted, paths, zrange):
        self._predicted = predicted
        self._paths = paths
        self._tar = tar
        self._zrange = zrange

    def __del__(self):
        self._tar.close()

    def paths(self, index):
        assert(index >= 0 and index < len(self))
        return iter(self._paths[index])

    def __getitem__(self, index):
        return self._load_block(index)

    def __len__(self):
        return len(self._paths)

    def __iter__(self):
        for i in range(len(self)):
            yield self._load_block(i)

    def _load_block(self, index):
        from dials.array_family import flex
        partials = [self._load_partial(p) for p in self.paths(index)]
        part = flex.partial_shoebox()
        ind = flex.int()
        for i, s in partials:
            part.extend(s)
            ind.extend(i)
        partials = (ind, part)
        indices, shoeboxes = zip(*self._merge_partials(partials))
        return (self._predicted.select(flex.size_t(indices)), shoeboxes)

    def _load_partial(self, path):
        import cPickle as pickle
        return pickle.load(self._tar.extractfile(path))

    def _merge_partials(self, partials):
        from collections import defaultdict
        from dials.array_family import flex
        tomerge = defaultdict(list)
        for index, shoebox in zip(*partials):
            tomerge[index].append(shoebox)

        return [(i, flex.partial_shoebox(m).merge(self._zrange)) for i, m in tomerge.iteritems()]


class ShoeboxSplitter(object):

    def __init__(self, filename, blocks, predicted):
        from collections import defaultdict
        import tarfile
        assert(len(blocks) > 1)
        assert(blocks[0] == 0)
        assert(all(b >= a for a, b in zip(blocks, blocks[1:])))
        self._predicted = predicted
        self._blocks = blocks
        self._tar = tarfile.open(filename, 'w')
        self._paths = defaultdict(list)
        self._lookup = []
        count = 0
        for i in range(max(blocks)):
            if i >= blocks[count+1]:
                count += 1
            self._lookup.append(count)

    def __del__(self):
        self._dump_index()
        self._tar.close()

    def __setitem__(self, index, item):
        self._dump_block(index, *item)

    def __len__(self):
        return len(self._blocks) - 1

    def block(self, frame):
        return self._lookup[frame]

    def _dump_block(self, index, indices, shoeboxes):
        from scitbx.array_family import flex
        assert(len(indices) == len(shoeboxes))
        for block, sind in self._split(indices).iteritems():
            sind = flex.size_t(sind)
            self._dump_shoeboxes(block, indices.select(sind),
                shoeboxes.select(sind))

    def _split(self, indices):
        from collections import defaultdict
        sperb = defaultdict(list)
        for i in range(len(indices)):
            sperb[self.block(int(self._predicted[indices[i]].frame_number))].append(i)
        return sperb

    def _dump_shoeboxes(self, block, indices, shoeboxes):
        import cPickle as pickle
        from StringIO import StringIO
        from uuid import uuid4
        from time import time
        data = StringIO(pickle.dumps((indices, shoeboxes),
            protocol=pickle.HIGHEST_PROTOCOL))
        info = self._tar.tarinfo()
        info.name = '%s.p' % uuid4()
        info.size = data.len
        info.mtime = int(time())
        self._tar.addfile(info, data)
        self._paths[block].append(info.name)

    def _dump_index(self):
        import cPickle as pickle
        from StringIO import StringIO
        from time import time
        zrange = min(self._blocks), max(self._blocks)
        index = { 'paths' : self._paths, 'predicted' : self._predicted, 'zrange' : zrange }
        data = StringIO(pickle.dumps(index, protocol=pickle.HIGHEST_PROTOCOL))
        info = self._tar.tarinfo()
        info.name = 'index.p'
        info.size = data.len
        info.mtime = int(time())
        self._tar.addfile(info, data)


def load_partial_shoeboxes(filename):
    import tarfile
    import cPickle as pickle
    tar = tarfile.open(filename)
    index = pickle.load(tar.extractfile('index.p'))
    return PartialShoeboxMerger(tar, index['predicted'], index['paths'], index['zrange'])

if __name__ == '__main__':

    import tarfile
    import StringIO
    from time import time
    import cPickle as pickle
    from scitbx.array_family import flex

    predicted = flex.int(100, 0)

    paths = [[] for i in range(10)]

    outfile = tarfile.open("temp.tar", 'w')
    for i in range(100):
        data = StringIO.StringIO(pickle.dumps((i // 10, i)))
        info = outfile.tarinfo()
        info.name = '%04x.p' % i
        info.size = data.len
        info.mtime = int(time())
        outfile.addfile(info, data)
        predicted[i] = i
        paths[i % 10].append(info.name)


    index = { 'paths' : paths, 'predicted' : predicted }
    data = StringIO.StringIO(pickle.dumps(index))
    info = outfile.tarinfo()
    info.name = 'index.p'
    info.size = data.len
    info.mtime = int(time())
    outfile.addfile(info, data)

    outfile.close()


    merger = load_partial_shoeboxes('temp.tar')

    for m in merger:
        p, s = m
        print list(p), s
