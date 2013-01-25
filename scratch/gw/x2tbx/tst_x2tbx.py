from __future__ import division
import scitbx.array_family.flex # explicit import
import cctbx.uctbx # explicit import

def tst_x2tbx(mtz_file):
    import x2tbx
    from iotbx import mtz

    mtz_obj = mtz.object(mtz_file)

    i_data = None
    sigi_data = None

    mi = mtz_obj.extract_miller_indices()

    unit_cell = None

    for crystal in mtz_obj.crystals():
        unit_cell = crystal.unit_cell()
        for dataset in crystal.datasets():
            for column in dataset.columns():
                if column.label() == 'I':
                    i_data = column.extract_values(
                        not_a_number_substitute = 0.0)
                if column.label() == 'SIGI':
                    sigi_data = column.extract_values(
                        not_a_number_substitute = 0.0)

    assert(i_data)
    assert(sigi_data)

    r = x2tbx.ReflectionList()
    r.setup(mi, i_data, sigi_data)
    r.set_unit_cell(unit_cell.parameters())
    r.merge()
    print r.i_sigma()
    print r.rmerge()

    indices = r.get_indices()
    print len(indices), len(mi)

    n_shells = 20

    r.setup_resolution_shells(n_shells)

    high = r.shell_high_limits()
    low = r.shell_low_limits()

    rmerges = r.rmerge_shells()
    isigmas = r.i_sigma_shells()
    tisigmas = r.total_i_sigma_shells()

    n_tot = 0

    for j in range(n_shells):
        shell = r.get_shell(j)
        print '%.3f %6.3f %4d %.3f %6.2f %6.2f' % (
            high[j], low[j], len(shell), rmerges[j], isigmas[j], tisigmas[j])
        n_tot += len(shell)

    assert(n_tot == len(indices))

    print 'OK'

if __name__ == '__main__':
    import sys
    tst_x2tbx(sys.argv[1])
