from __future__ import annotations

import socket
import subprocess
import sys
import time
import timeit
import urllib.request
from xml.dom import minidom

import procrunner
import pytest


@pytest.fixture
def server_port(tmp_path) -> int:
    """Fixture to load a find_spots_server server"""
    if sys.hexversion >= 0x3080000 and sys.platform == "darwin":
        pytest.skip("find_spots server known to be broken on MacOS with Python 3.8+")
    if sys.platform == "win32":
        pytest.skip("find_spots server is not supported on Windows")

    # Find a free port to run the server on
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("", 0))
        host, port = sock.getsockname()
    # Start the server
    server_command = ["dials.find_spots_server", f"port={port}", "nproc=3"]
    p = subprocess.Popen(server_command, cwd=tmp_path)
    wait_for_server(port)
    yield port
    p.terminate()
    p.wait(timeout=3)
    p.kill()


def test_server_return_codes(dials_data, server_port):
    try:
        first_file = sorted(
            dials_data("centroid_test_data", pathlib=True).glob("*.cbf")
        )[0]
        response = urllib.request.urlopen(
            f"http://127.0.0.1:{server_port}/{first_file}"
        )
        assert response.code == 200
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(f"http://127.0.0.1:{server_port}/some/junk/filename")
    finally:
        result = procrunner.run(
            ["dials.find_spots_client", f"port={server_port}", "stop"]
        )
        assert not result.returncode and not result.stderr


def test_find_spots_server_client(dials_data, tmp_path, server_port):
    filenames = sorted(dials_data("centroid_test_data", pathlib=True).glob("*.cbf"))
    try:
        exercise_client(server_port=server_port, filenames=filenames)
    finally:
        result = procrunner.run(
            ["dials.find_spots_client", f"port={server_port}", "stop"]
        )
        assert not result.returncode and not result.stderr


def wait_for_server(port, max_wait=20):
    print(f"Waiting up to {max_wait} seconds for server to start")
    server_ok = False
    start_time = timeit.default_timer()
    max_time = start_time + max_wait
    while (timeit.default_timer() < max_time) and not server_ok:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.connect(("127.0.0.1", port))
            server_ok = True
        except OSError as e:
            if (e.errno != 111) and (e.errno != 61):
                raise
            # ignore connection failures (111 connection refused on linux; 61 connection refused on mac)
            time.sleep(0.1)
    if not server_ok:
        raise Exception("Server failed to start after %d seconds" % max_wait)
    print(
        "dials.find_spots_server up after %f seconds"
        % (timeit.default_timer() - start_time)
    )


def exercise_client(server_port, filenames):
    assert filenames
    client_command = [
        "dials.find_spots_client",
        f"port={server_port}",
        "min_spot_size=3",
        "algorithm=dispersion",
        "nproc=1",
        filenames[0],
    ]

    index_client_command = client_command + [
        "index=True",
        "indexing.method=fft1d",
        "max_refine=10",
    ]
    print(index_client_command)
    result = procrunner.run(index_client_command)
    assert not result.returncode and not result.stderr
    out = f"<document>{result.stdout}</document>"

    xmldoc = minidom.parseString(out)
    assert len(xmldoc.getElementsByTagName("image")) == 1
    assert len(xmldoc.getElementsByTagName("spot_count")) == 1
    assert len(xmldoc.getElementsByTagName("spot_count_no_ice")) == 1
    assert len(xmldoc.getElementsByTagName("d_min")) == 1
    assert len(xmldoc.getElementsByTagName("total_intensity")) == 1
    assert len(xmldoc.getElementsByTagName("unit_cell")) == 1
    assert len(xmldoc.getElementsByTagName("n_indexed")) == 1
    assert len(xmldoc.getElementsByTagName("fraction_indexed")) == 1

    unit_cell = [
        float(f)
        for f in xmldoc.getElementsByTagName("unit_cell")[0].childNodes[0].data.split()
    ]

    assert unit_cell == pytest.approx(
        [39.90, 42.67, 42.37, 89.89, 90.10, 90.13], abs=1e-1
    )

    client_command = client_command + filenames[1:]
    result = procrunner.run(client_command)
    assert not result.returncode and not result.stderr
    out = f"<document>{result.stdout}</document>"

    xmldoc = minidom.parseString(out)
    images = xmldoc.getElementsByTagName("image")
    assert len(images) == 9
    spot_counts = sorted(
        int(node.childNodes[0].data)
        for node in xmldoc.getElementsByTagName("spot_count")
    )
    assert spot_counts == sorted([203, 196, 205, 209, 195, 205, 203, 207, 189])
    spot_counts_no_ice = sorted(
        int(node.childNodes[0].data)
        for node in xmldoc.getElementsByTagName("spot_count_no_ice")
    )
    assert spot_counts_no_ice == sorted([169, 171, 175, 176, 177, 184, 193, 195, 196])
    d_min = sorted(
        float(node.childNodes[0].data) for node in xmldoc.getElementsByTagName("d_min")
    )
    assert d_min == sorted([1.45, 1.47, 1.55, 1.55, 1.56, 1.59, 1.61, 1.61, 1.64])
