from __future__ import division

def resolutionizer(mtz_file, n_shells = 30):
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

    r.setup_resolution_shells(n_shells)

    high = r.shell_high_limits()
    low = r.shell_low_limits()

    rmerges = r.rmerge_shells()
    isigmas = r.i_sigma_shells()
    tisigmas = r.total_i_sigma_shells()

    print 'High   Low    Nref  Rmerge Mn(I/s) I/s'
    for j in range(n_shells):
        shell = r.get_shell(j)
        print '%.3f %6.3f %5d %.3f %6.2f %6.2f' % (
            high[j], low[j], len(shell), rmerges[j], isigmas[j], tisigmas[j])

    return

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        resolutionizer(sys.argv[1])
    else:
        resolutionizer(sys.argv[1], n_shells = int(sys.argv[2]))
