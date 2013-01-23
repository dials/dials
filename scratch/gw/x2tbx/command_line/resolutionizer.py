from __future__ import division
import scitbx.array_family.flex # explicit import
import cctbx.uctbx # explicit import

def resolutionizer(mtz_file):
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

    r = x2tbx.resolutionizer()
    r.set_unit_cell(unit_cell.parameters())
    r.setup(mi, i_data, sigi_data)
    print r.isig()

if __name__ == '__main__':
    import sys
    resolutionizer(sys.argv[1])
