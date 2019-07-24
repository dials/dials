from __future__ import absolute_import, division, print_function

from future import standard_library

standard_library.install_aliases()
import httplib
import json
import os
import socket as pysocket
import sys
import urllib.error
import urllib.parse
import urllib.request

import libtbx.phil


def work(host, port, filename, params):
    conn = httplib.HTTPConnection(host, port)
    path = filename
    for param in params:
        path += ";%s" % param
    conn.request("GET", path)
    return conn.getresponse().read()


def _nproc():
    from libtbx.introspection import number_of_processors

    return number_of_processors(return_value_if_unknown=-1)


def response_to_xml(d):

    if "n_spots_total" in d:
        response = (
            """\
<image>%(image)s</image>
<spot_count>%(n_spots_total)s</spot_count>
<spot_count_no_ice>%(n_spots_no_ice)s</spot_count_no_ice>
<d_min>%(estimated_d_min).2f</d_min>
<d_min_method_1>%(d_min_distl_method_1).2f</d_min_method_1>
<d_min_method_2>%(d_min_distl_method_2).2f</d_min_method_2>
<total_intensity>%(total_intensity).0f</total_intensity>"""
            % d
        )

    else:
        assert "error" in d
        return "<response>\n%s\n</response>" % d["error"]

    if "lattices" in d:
        from dxtbx.model.crystal import CrystalFactory

        for lattice in d["lattices"]:
            crystal = CrystalFactory.from_dict(lattice["crystal"])
            response = "\n".join(
                [
                    response,
                    "<unit_cell>%.6g %.6g %.6g %.6g %.6g %.6g</unit_cell>"
                    % (crystal.get_unit_cell().parameters()),
                ]
            )
        response = "\n".join(
            [
                response,
                "<n_indexed>%i</n_indexed>" % d["n_indexed"],
                "<fraction_indexed>%.2f</fraction_indexed>" % d["fraction_indexed"],
            ]
        )
    if "integrated_intensity" in d:
        response = "\n".join(
            [
                response,
                "<integrated_intensity>%.0f</integrated_intensity>"
                % d["integrated_intensity"],
            ]
        )

    return "<response>\n%s\n</response>" % response


def work_all(
    host,
    port,
    filenames,
    params,
    plot=False,
    table=False,
    json_file=None,
    grid=None,
    nproc=None,
):
    from multiprocessing.pool import ThreadPool as thread_pool

    if nproc is None:
        nproc = _nproc()
    pool = thread_pool(processes=nproc)
    threads = {}
    for filename in filenames:
        threads[filename] = pool.apply_async(work, (host, port, filename, params))
    results = []
    for filename in filenames:
        response = threads[filename].get()
        d = json.loads(response)
        results.append(d)
        print(response_to_xml(d))

    if json_file is not None:
        "Writing results to %s" % json_file
        with open(json_file, "wb") as f:
            json.dump(results, f)

    if plot or table:
        from scitbx.array_family import flex
        from libtbx import group_args
        from dials.algorithms.spot_finding.per_image_analysis import (
            plot_stats,
            print_table,
        )

        estimated_d_min = flex.double()
        d_min_distl_method_1 = flex.double()
        d_min_distl_method_2 = flex.double()
        n_spots_total = flex.int()
        n_spots_no_ice = flex.int()
        total_intensity = flex.double()

        for d in results:
            estimated_d_min.append(d["estimated_d_min"])
            d_min_distl_method_1.append(d["d_min_distl_method_1"])
            d_min_distl_method_2.append(d["d_min_distl_method_2"])
            n_spots_total.append(d["n_spots_total"])
            n_spots_no_ice.append(d["n_spots_no_ice"])
            total_intensity.append(d["total_intensity"])

        stats = group_args(
            n_spots_total=n_spots_total,
            n_spots_no_ice=n_spots_no_ice,
            n_spots_4A=None,
            total_intensity=total_intensity,
            estimated_d_min=estimated_d_min,
            d_min_distl_method_1=d_min_distl_method_1,
            d_min_distl_method_2=d_min_distl_method_2,
            noisiness_method_1=None,
            noisiness_method_2=None,
        )

        if plot:
            plot_stats(stats)
        if table:
            print_table(stats)

        if grid is not None:
            from matplotlib import pyplot

            n_spots_no_ice.reshape(flex.grid(grid))
            print(n_spots_no_ice.size())

            pyplot.figure()
            pyplot.pcolormesh(n_spots_no_ice.as_numpy_array(), cmap=pyplot.cm.Reds)
            pyplot.savefig("spot_count.png")


def stop(host, port, nproc):
    stopped = 0
    for j in range(nproc):
        try:
            url_request = urllib.request.Request("http://%s:%s/Ctrl-C" % (host, port))
            socket = urllib.request.urlopen(url_request, None, 3)
            if socket.getcode() == "200":
                stopped = stopped + 1
            else:
                print("socket returned code", socket.getcode())
        except (pysocket.timeout, urllib.error.HTTPError) as e:
            print("error on stopping server:", e)
        except urllib.error.URLError as e:
            if e.reason.errno != 111:
                print("error on stopping server:", e)
        except pysocket.error:
            # Assuming this means the server killed itself before the reply left the send buffer.
            stopped = stopped + 1
        except httplib.BadStatusLine:
            # Regular occurrence. Probably means the server stopped anyway.
            stopped = stopped + 1
    return stopped


phil_scope = libtbx.phil.parse(
    """\
nproc = Auto
  .type = int(value_min=1)
host = localhost
  .type = str
port = 1701
  .type = int(value_min=1)
plot = False
  .type = bool
table = False
  .type = bool
json = None
  .type = path
grid = None
  .type = ints(size=2, value_min=1)
"""
)

if __name__ == "__main__":
    import select

    args = sys.argv[1:]

    if os.name != "nt":
        r, w, x = select.select([sys.stdin], [], [], 0)
        if len(r) > 0:
            args.extend([l.strip() for rr in r for l in rr.readlines()])

    filenames = [arg for arg in args if os.path.isfile(arg)]
    args = [arg for arg in args if arg not in filenames]

    interp = phil_scope.command_line_argument_interpreter()
    params, unhandled = interp.process_and_fetch(
        args, custom_processor="collect_remaining"
    )
    params = params.extract()

    if params.nproc is libtbx.Auto:
        nproc = None
        params.nproc = 1024
    else:
        nproc = params.nproc

    if len(unhandled) and unhandled[0] == "stop":
        stopped = stop(params.host, params.port, params.nproc)
        print("Stopped %d findspots processes" % stopped)
    elif len(unhandled) and unhandled[0] == "ping":
        url = "http://%s:%i" % (params.host, params.port)
        try:
            data = urllib.request.urlopen(url).read()
            print("Success")
            sys.exit(0)
        except Exception:
            print("Failure")
            sys.exit(1)
    else:
        if len(filenames) == 1:
            response = work(params.host, params.port, filenames[0], unhandled)

            print(response_to_xml(json.loads(response)))
        else:
            work_all(
                params.host,
                params.port,
                filenames,
                unhandled,
                plot=params.plot,
                table=params.table,
                json_file=params.json,
                grid=params.grid,
                nproc=nproc,
            )
