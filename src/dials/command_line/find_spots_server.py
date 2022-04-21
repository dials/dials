from __future__ import annotations

import http.server as server_base
import json
import logging
import multiprocessing
import sys
import time
import urllib.parse

import libtbx.phil
from cctbx import uctbx
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.introspection import number_of_processors

from dials.algorithms.indexing import indexer
from dials.algorithms.integration.integrator import create_integrator
from dials.algorithms.profile_model.factory import ProfileModelFactory
from dials.algorithms.spot_finding import per_image_analysis
from dials.array_family import flex
from dials.command_line.find_spots import phil_scope as find_spots_phil_scope
from dials.command_line.index import phil_scope as index_phil_scope
from dials.command_line.integrate import phil_scope as integrate_phil_scope
from dials.util import Sorry, show_mail_handle_errors
from dials.util.options import ArgumentParser

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


def _filter_by_resolution(experiments, reflections, d_min=None, d_max=None):
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    d_star_sq = flex.pow2(reflections["rlp"].norms())
    reflections["d"] = uctbx.d_star_sq_as_d(d_star_sq)
    # Filter based on resolution
    if d_min is not None:
        selection = reflections["d"] >= d_min
        reflections = reflections.select(selection)
        logger.debug(f"Selected {len(reflections)} reflections with d >= {d_min:f}")

    # Filter based on resolution
    if d_max is not None:
        selection = reflections["d"] <= d_max
        reflections = reflections.select(selection)
        logger.debug(f"Selected {len(reflections)} reflections with d <= {d_max:f}")
    return reflections


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
    interp = phil_scope.command_line_argument_interpreter()
    params, unhandled = interp.process_and_fetch(
        cl, custom_processor="collect_remaining"
    )
    filter_ice = params.extract().ice_rings.filter
    ice_rings_width = params.extract().ice_rings.width
    index = params.extract().index
    integrate = params.extract().integrate
    indexing_min_spots = params.extract().indexing_min_spots

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
    if params.spotfinder.scan_range and len(experiments) > 1:
        # This means we've imported a sequence of still image: select
        # only the experiment, i.e. image, we're interested in
        ((start, end),) = params.spotfinder.scan_range
        experiments = experiments[start - 1 : end]

    # Avoid overhead of calculating per-pixel resolution masks in spotfinding
    # and instead perform post-filtering of spot centroids by resolution
    d_min = params.spotfinder.filter.d_min
    d_max = params.spotfinder.filter.d_max
    params.spotfinder.filter.d_min = None
    params.spotfinder.filter.d_max = None

    t0 = time.perf_counter()
    reflections = flex.reflection_table.from_observations(experiments, params)

    if d_min or d_max:
        reflections = _filter_by_resolution(
            experiments, reflections, d_min=d_min, d_max=d_max
        )

    t1 = time.perf_counter()
    logger.info("Spotfinding took %.2f seconds", t1 - t0)

    imageset = experiments.imagesets()[0]
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    stats = per_image_analysis.stats_for_reflection_table(
        reflections, filter_ice=filter_ice, ice_rings_width=ice_rings_width
    )._asdict()
    t2 = time.perf_counter()
    logger.info("Resolution analysis took %.2f seconds", t2 - t1)

    if index and stats["n_spots_no_ice"] > indexing_min_spots:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)

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
            and imageset.get_scan().is_still()
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
            t3 = time.perf_counter()
            logger.info("Indexing took %.2f seconds", t3 - t2)

        if integrate and "lattices" in stats:
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
                        "Warning: %d reference spots were not matched to predictions",
                        len(reference) - matched.count(True),
                    )
                    logger.info("*" * 80)
                    logger.info("")

                # Compute the profile model
                experiments = ProfileModelFactory.create(params, experiments, reference)

                # Compute the bounding box
                predicted.compute_bbox(experiments)

                # Create the integrator
                integrator = create_integrator(params, experiments, predicted)

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
                t4 = time.perf_counter()
                logger.info("Integration took %.2f seconds", t4 - t3)

    return stats


class handler(server_base.BaseHTTPRequestHandler):
    def do_GET(self):
        """Respond to a GET request."""
        if self.path == "/Ctrl-C":
            self.send_response(200)
            self.end_headers()

            global stop
            stop = True
            return

        filename = self.path.split(";")[0]
        params = self.path.split(";")[1:]

        # If we're passing a url through, then unquote and ignore leading /
        if "%3A//" in filename:
            filename = urllib.parse.unquote(filename[1:])

        d = {"image": filename}

        try:
            stats = work(filename, params)
            d.update(stats)
            response = 200
        except Exception as e:
            d["error"] = str(e)
            response = 500

        self.send_response(response)
        self.send_header("Content-type", "application/json")
        self.end_headers()
        response = json.dumps(d).encode()
        self.wfile.write(response)


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
        proc = multiprocessing.Process(target=serve, args=(httpd,))
        proc.daemon = True
        proc.start()
    serve(httpd)
    httpd.server_close()
    print(time.asctime(), "done")


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.find_spots_server [options]"

    # Python 3.8 on macOS... needs fork
    if sys.hexversion >= 0x3080000 and sys.platform == "darwin":
        multiprocessing.set_start_method("fork")

    parser = ArgumentParser(usage=usage, phil=phil_scope, epilog=help_message)
    params, options = parser.parse_args(args, show_diff_phil=True)
    if params.nproc is libtbx.Auto:
        params.nproc = number_of_processors(return_value_if_unknown=-1)
    main(params.nproc, params.port)


if __name__ == "__main__":
    run()
