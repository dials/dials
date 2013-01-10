from toytbx import make_list

def tst_toytbx():
    assert(make_list(4) == [j for j in range(4)])
    print 'OK'

if __name__ == '__main__':
    tst_toytbx()


