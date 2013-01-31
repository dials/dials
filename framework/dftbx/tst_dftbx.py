def tst_dftbx():
    import dftbx
    int_ol = dftbx.int_ol()
    for j in range(1000):
        int_ol.push_back(j)
    assert(int_ol.size() == 1000)
    print 'OK'

if __name__ == '__main__':
    tst_dftbx()



