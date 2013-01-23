import scitbx.array_family.flex
import x2tbx

def tst_x2tbx():
    assert(x2tbx.make_list(4) == [j for j in range(4)])
    assert(sum(x2tbx.make_flex(10)) == x2tbx.sum(x2tbx.make_flex(10)))
    print 'OK'

if __name__ == '__main__':
    tst_x2tbx()


