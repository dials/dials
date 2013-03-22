import scitbx.array_family.flex
import toytbx

def tst_toytbx():
    assert(toytbx.make_list(4) == [j for j in range(4)])
    assert(sum(toytbx.make_flex(10)) == toytbx.sum(toytbx.make_flex(10)))
    print 'OK'

if __name__ == '__main__':
    tst_toytbx()
