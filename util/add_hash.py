from __future__ import absolute_import, division, print_function


def enhash(e, h, k, l):
    return e * (2 ** 30) + (h + 512) * (2 ** 20) + (k + 512) * (2 ** 10) + (l + 512)


def dehash(hash_value):
    if isinstance(hash_value, type(42)):
        e = hash_value // 2 ** 30
        h = (hash_value % (2 ** 30)) // (2 ** 20) - 512
        k = (hash_value % (2 ** 20)) // (2 ** 10) - 512
        l = (hash_value % (2 ** 10)) - 512
    else:
        e = hash_value / 2 ** 30
        h = (hash_value % (2 ** 30)) / (2 ** 20) - 512
        k = (hash_value % (2 ** 20)) / (2 ** 10) - 512
        l = (hash_value % (2 ** 10)) - 512
    return e, h, k, l


def add_hash(integrated_data):
    """Add hash = 2^30 * entering + 2^20 * (h+512) + 2^10 * (k+512) + (l+512)
    as new column to reflection table - should be P1 unique for 360 degree
    scans"""

    from dials.array_family import flex

    integrated_data = integrated_data.select(integrated_data["id"] >= 0)
    assert max(integrated_data["id"]) == 0

    h, k, l = integrated_data["miller_index"].as_vec3_double().parts()
    h = h.iround()
    k = k.iround()
    l = l.iround()
    e = integrated_data["entering"].as_int()

    hash_value = enhash(e, h, k, l)
    integrated_data["hash"] = hash_value

    return integrated_data


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise RuntimeError("%s strong.refl" % sys.argv[0])

    from dials.array_family import flex

    integrated_data = flex.reflection_table.from_file(sys.argv[1])

    # keep only flagged as integrated reflections
    sel = integrated_data.get_flags(integrated_data.flags.integrated)
    integrated_data = integrated_data.select(sel)

    # only fully recorded reflections
    sel = integrated_data["partiality"] > 0.99
    integrated_data = integrated_data.select(sel)

    integrated_data = add_hash(integrated_data)
    hash_value = integrated_data["hash"]

    for h in hash_value:
        sel = hash_value == h
        assert sel.count(True) == 1

    print(flex.min(hash_value), flex.max(hash_value))
