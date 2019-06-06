from __future__ import absolute_import, division, print_function

import BaseHTTPServer as server_base
import logging
import os
import sys
import time
from multiprocessing import Process

import libtbx.load_env
import libtbx.phil

from dials.util import Sorry

logger = logging.getLogger("dials.command_line.find_spots_server")

help_message = """\
A client/server version of dials.find_spots with additional analysis including
estimation of resolution limits. Intended for quick feedback of image quality
during grid scans and data collections.

On the server machine::

  dials.find_spots_server [nproc=8] [port=1234]

On the client machine::

  dials.find_spots_client [host=hostname] [port=1234] [nproc=8] /path/to/image.cbf

The client will return a short xml string indicating the number of spots found
and several estimates of the resolution limit.

e.g.::

  <response>
  <image>/path/to/image_0001.cbf</image>
  <spot_count>352</spot_count>
  <spot_count_no_ice>263</spot_count_no_ice>
  <d_min>1.46</d_min>
  <d_min_method_1>1.92</d_min_method_1>
  <d_min_method_2>1.68</d_min_method_2>
  <total_intensity>56215</total_intensity>
  </response>


* ``spot_count`` is the total number of spots found in given image
* ``spot_count_no_ice`` is the number of spots found excluding those at resolutions
  where ice rings may be found
* ``d_min_method_1`` is equivalent to distl's resolution estimate method 1
* ``d_min_method_2`` is equivalent to distl's resolution estimate method 2
* ``total_intensity`` is the total intensity of all strong spots excluding those
  at resolutions where ice rings may be found

Any valid ``dials.find_spots`` parameter may be passed to
``dials.find_spots_client``, e.g.::

  dials.find_spots_client /path/to/image.cbf min_spot_size=2 d_min=2

To stop the server::

  dials.find_spots_client stop [host=hostname] [port=1234]

"""

stop = False


def work(filename, cl=None):
    if cl is None:
        cl = []

    phil_scope = libtbx.phil.parse(
        """\
ice_rings {
  filter = True
    .type = bool
  width = 0.004
    .type = float(value_min=0.0)
}
index = False
  .type = bool
integrate = False
  .type = bool
indexing_min_spots = 10
  .type = int(value_min=1)
"""
    )
    if not os.access(filename, os.R_OK):
        raise RuntimeError("Server does not have read access to file %s" % filename)
    interp = phil_scope.command_line_argument_interpreter()
    params, unhandled = interp.process_and_fetch(
        cl, custom_processor="collect_remaining"
    )
    filter_ice = params.extract().ice_rings.filter
    ice_rings_width = params.extract().ice_rings.width
    index = params.extract().index
    integrate = params.extract().integrate
    indexing_min_spots = params.extract().indexing_min_spots

    from dials.command_line.find_spots import phil_scope as find_spots_phil_scope
    from dxtbx.model.experiment_list import ExperimentListFactory
    from dials.array_family import flex

    interp = find_spots_phil_scope.command_line_argument_interpreter()
    phil_scope, unhandled = interp.process_and_fetch(
        unhandled, custom_processor="collect_remaining"
    )
    logger.info("The following spotfinding parameters have been modified:")
    logger.info(find_spots_phil_scope.fetch_diff(source=phil_scope).as_str())
    params = phil_scope.extract()
    # no need to write the hot mask in the server/client
    params.spotfinder.write_hot_mask = False
    experiments = ExperimentListFactory.from_filenames([filename])
    t0 = time.time()
    reflections = flex.reflection_table.from_observations(experiments, params)
    t1 = time.time()
    logger.info("Spotfinding took %.2f seconds" % (t1 - t0))
    from dials.algorithms.spot_finding import per_image_analysis

    imageset = experiments.imagesets()[0]
    scan = imageset.get_scan()
    if scan is not None:
        i = scan.get_array_range()[0]
    else:
        i = 0
    stats = per_image_analysis.stats_single_image(
        imageset,
        reflections,
        i=i,
        plot=False,
        filter_ice=filter_ice,
        ice_rings_width=ice_rings_width,
    )
    stats = stats.__dict__
    t2 = time.time()
    logger.info("Resolution analysis took %.2f seconds" % (t2 - t1))

    if index and stats["n_spots_no_ice"] > indexing_min_spots:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)
        from dials.algorithms.indexing import indexer
        from dials.command_line.index import phil_scope as index_phil_scope

        interp = index_phil_scope.command_line_argument_interpreter()
        phil_scope, unhandled = interp.process_and_fetch(
            unhandled, custom_processor="collect_remaining"
        )
        logger.info("The following indexing parameters have been modified:")
        index_phil_scope.fetch_diff(source=phil_scope).show()
        params = phil_scope.extract()

        if (
            imageset.get_goniometer() is not None
            and imageset.get_scan() is not None
            and imageset.get_scan().get_oscillation()[1] == 0
        ):
            imageset.set_goniometer(None)
            imageset.set_scan(None)

        try:
            idxr = indexer.Indexer.from_parameters(
                reflections, experiments, params=params
            )
            indexing_results = []
            idxr.index()
            indexed_sel = idxr.refined_reflections.get_flags(
                idxr.refined_reflections.flags.indexed
            )
            indexed_sel &= ~(
                idxr.refined_reflections.get_flags(
                    idxr.refined_reflections.flags.centroid_outlier
                )
            )
            for i_expt, expt in enumerate(idxr.refined_experiments):
                sel = idxr.refined_reflections["id"] == i_expt
                sel &= indexed_sel
                indexing_results.append(
                    {
                        "crystal": expt.crystal.to_dict(),
                        "n_indexed": sel.count(True),
                        "fraction_indexed": sel.count(True) / sel.size(),
                    }
                )
            stats["lattices"] = indexing_results
            stats["n_indexed"] = indexed_sel.count(True)
            stats["fraction_indexed"] = indexed_sel.count(True) / len(reflections)
        except Exception as e:
            logger.error(e)
            stats["error"] = str(e)
        finally:
            t3 = time.time()
            logger.info("Indexing took %.2f seconds" % (t3 - t2))

        if integrate and "lattices" in stats:

            from dials.algorithms.profile_model.factory import ProfileModelFactory
            from dials.algorithms.integration.integrator import IntegratorFactory
            from dials.command_line.integrate import phil_scope as integrate_phil_scope

            interp = integrate_phil_scope.command_line_argument_interpreter()
            phil_scope, unhandled = interp.process_and_fetch(
                unhandled, custom_processor="collect_remaining"
            )
            logger.error("The following integration parameters have been modified:")
            integrate_phil_scope.fetch_diff(source=phil_scope).show()
            params = phil_scope.extract()

            try:
                params.profile.gaussian_rs.min_spots = 0

                experiments = idxr.refined_experiments
                reference = idxr.refined_reflections

                predicted = flex.reflection_table.from_predictions_multi(
                    experiments,
                    dmin=params.prediction.d_min,
                    dmax=params.prediction.d_max,
                    margin=params.prediction.margin,
                    force_static=params.prediction.force_static,
                )

                matched, reference, unmatched = predicted.match_with_reference(
                    reference
                )
                assert len(matched) == len(predicted)
                assert matched.count(True) <= len(reference)
                if matched.count(True) == 0:
                    raise Sorry(
                        """
            Invalid input for reference reflections.
            Zero reference spots were matched to predictions
          """
                    )
                elif matched.count(True) != len(reference):
                    logger.info("")
                    logger.info("*" * 80)
                    logger.info(
                        "Warning: %d reference spots were not matched to predictions"
                        % (len(reference) - matched.count(True))
                    )
                    logger.info("*" * 80)
                    logger.info("")

                # Compute the profile model
                experiments = ProfileModelFactory.create(params, experiments, reference)

                # Compute the bounding box
                predicted.compute_bbox(experiments)

                # Create the integrator
                integrator = IntegratorFactory.create(params, experiments, predicted)

                # Integrate the reflections
                reflections = integrator.integrate()

                # print len(reflections)

                stats["integrated_intensity"] = flex.sum(
                    reflections["intensity.sum.value"]
                )
            except Exception as e:
                logger.error(e)
                stats["error"] = str(e)
            finally:
                t4 = time.time()
                logger.info("Integration took %.2f seconds" % (t4 - t3))

    return stats


class handler(server_base.BaseHTTPRequestHandler):
    def do_GET(s):
        """Respond to a GET request."""
        s.send_response(200)
        s.send_header("Content-type", "text/xml")
        s.end_headers()
        if s.path == "/Ctrl-C":
            global stop
            stop = True
            return
        filename = s.path.split(";")[0]
        params = s.path.split(";")[1:]

        d = {"image": filename}

        try:
            stats = work(filename, params)
            d.update(stats)

        except Exception as e:
            d["error"] = str(e)

        import json

        response = json.dumps(d)
        s.wfile.write(response)


def serve(httpd):
    try:
        while not stop:
            httpd.handle_request()
    except KeyboardInterrupt:
        pass


phil_scope = libtbx.phil.parse(
    """\
nproc = Auto
  .type = int(value_min=1)
port = 1701
  .type = int(value_min=1)
"""
)


def main(nproc, port):
    server_class = server_base.HTTPServer
    httpd = server_class(("", port), handler)
    print(time.asctime(), "Serving %d processes on port %d" % (nproc, port))

    for j in range(nproc - 1):
        proc = Process(target=serve, args=(httpd,))
        proc.daemon = True
        proc.start()
    serve(httpd)
    httpd.server_close()
    print(time.asctime(), "done")


if __name__ == "__main__":
    usage = "%s [options]" % libtbx.env.dispatcher_name

    from dials.util.options import OptionParser

    parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)
    params, options = parser.parse_args(show_diff_phil=True)
    if params.nproc is libtbx.Auto:
        from libtbx.introspection import number_of_processors

        params.nproc = number_of_processors(return_value_if_unknown=-1)
    main(params.nproc, params.port)
