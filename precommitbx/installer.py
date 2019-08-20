from __future__ import absolute_import, division, print_function

import contextlib
import multiprocessing.pool
import os
import stat
import sys
import time

import dials.precommitbx._precommitbx
import libtbx.introspection
import libtbx.load_env
import procrunner
import py
from future.moves.urllib.request import urlopen, Request
from tqdm import tqdm, trange

BOLD = "\033[1m"
GREEN = "\033[32m"
MAGENTA = "\033[1;35m"
NC = "\033[0m"
RED = "\033[1;31m"
YELLOW = "\033[1;33m"

precommit_home = py.path.local(abs(libtbx.env.build_path)).join("precommitbx")
precommitbx_version = dials.precommitbx._precommitbx.__version__

libffi_source_size = 940837
libffi_source_version = "3.2.1"
openssl_API_version = "(1, 0, 2, 19, 15)"  # import ssl; print(ssl._OPENSSL_API_VERSION)
openssl_source_size = 5349149
openssl_source_version = "1.0.2s"
python_source_size = 23017663
python_source_version = "3.7.4"
sqlite_source_size = 2833613
sqlite_source_version = 3290000
sqlite_version = "3.29.0"

environment_override = {
    "LD_LIBRARY_PATH": "",
    "PRE_COMMIT_HOME": precommit_home.join("cache"),
    "PYTHONPATH": "",
}

repo_prefix = "  {:.<15}:"
repo_no_precommit = "(no pre-commit hooks)"
repo_precommit_installed = GREEN + "pre-commit installed" + NC
repo_precommit_conflict = (
    RED + "pre-commit available but a different pre-commit hook is installed" + NC
)
repo_precommit_legacy = YELLOW + "pre-commit out of date" + NC
repo_precommit_available = MAGENTA + "pre-commit available but not installed" + NC


def clean_run(*args, **kwargs):
    stop_on_error = kwargs.pop("stop_on_error", False)
    env_override = environment_override.copy()
    env_override.update(kwargs.pop("environment_override", {}))
    result = procrunner.run(
        *args,
        environment_override=env_override,
        print_stderr=False,
        print_stdout=False,
        **kwargs
    )
    if stop_on_error and result.returncode:
        if result.stdout:
            print("\n".join(result.stdout.split("\n")[-10:]))
            print("---")
        if result.stderr:
            print("\n".join(result.stderr.split("\n")[-10:]))
            print("---")
        sys.exit(stop_on_error)
    return result


class ProgressOverall(object):
    def __init__(self):
        self._status = {}
        self._bar = None
        for package, steps in {
            "Python": 5,
            "OpenSSL": 5,
            "libffi": 5,
            "SQLite": 5,
            "Precommitbx": 1,
        }.items():
            self._status[package] = {"total": steps, "complete": 0, "partial": 0}

    def start(self):
        if self._bar:
            return
        self._bar = tqdm(
            bar_format="{l_bar}{bar}",
            desc="Setting up precommitbx",
            leave=True,
            smoothing=0.1,
            total=100 * sum(x["total"] for x in self._status.values()),
        )

    def done_part(self, package, partial):
        if self._status[package]["complete"] >= self._status[package]["total"]:
            return
        self._status[package]["partial"] = partial
        self._update()

    def done_step(self, package):
        self._status[package]["partial"] = 0
        if self._status[package]["complete"] >= self._status[package]["total"]:
            precommit_home.join(".tqdm-record").write(
                "{package}: package already completed - skipping step".format(
                    package=package
                ),
                "a",
            )
            return
        self._status[package]["complete"] += 1
        self._update()

    def done_package(self, package):
        self._status[package]["partial"] = 0
        self._status[package]["complete"] = self._status[package]["total"]
        self._update()

    def _update(self):
        n = int(100 * sum(p["complete"] + p["partial"] for p in self._status.values()))
        if n > self._bar.n:
            self._bar.n = n
            self._bar.refresh()

    def done(self):
        self._bar.close()


ProgressOverall = ProgressOverall()


class Progress(object):
    def __init__(self, description, length, step=None):
        self.pbar = tqdm(
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{remaining}]",
            desc="{:<18}".format(description),
            leave=False,
            smoothing=0.1,
            total=length,
            unit="",
        )
        self.step = step

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.pbar.n != self.pbar.total:
            precommit_home.join(".tqdm-record").write(
                "{pbar.desc}: {pbar.n} out of {pbar.total}\n".format(pbar=self.pbar),
                "a",
            )
        self.pbar.close()
        if self.step:
            ProgressOverall.done_step(self.step)

    def increment(self, *args, **kwargs):
        self.pbar.update(1)
        if self.step:
            ProgressOverall.done_part(self.step, self.pbar.n / self.pbar.total)


def stop_with_error(message, result):
    if result.stdout:
        print("\n".join(result.stdout.split("\n")[-10:]))
    if result.stderr:
        print("\n".join(result.stderr.split("\n")[-10:]))
    sys.exit(message)


def precommitbx_template(python):
    return "\n".join(
        [
            "#!/bin/bash",
            "# File generated by precommitbx",
            "export LD_LIBRARY_PATH=",
            "export PYTHONPATH=",
            "export PRE_COMMIT_HOME=" + precommit_home.join("cache").strpath,
            "export PATH=" + python.dirname + os.pathsep + "$PATH",
            "if [ ! -f .pre-commit-config.yaml ]; then",
            "  echo No pre-commit configuration. Skipping pre-commit checks.",
            "  exit 0",
            "fi",
            r'if grep -q "language_version.\+libtbx.python" .pre-commit-config.yaml; then',
            "  echo Repository contains legacy pre-commit configuration. Skipping pre-commit checks.",
            "  exit 0",
            "fi",
            "if [ ! -e " + python.dirname + os.path.sep + python.basename + " ]; then",
            "  echo Precommitbx installation not found. Run libtbx.precommit to fix.",
            "  exit 0",
            "fi",
            python.basename + " -m _precommitbx.main",
        ]
    )


def install_precommitbx_hook(path, python):
    with path.join(".git", "hooks", "pre-commit").open("w") as fh:
        fh.write(precommitbx_template(python))
        mode = os.fstat(fh.fileno()).st_mode
        mode |= stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        os.fchmod(fh.fileno(), stat.S_IMODE(mode))


def check_precommitbx_hook(path, python):
    hookfile = path.join(".git", "hooks", "pre-commit")
    if not hookfile.check():
        return False
    if not os.access(hookfile.strpath, os.X_OK):
        return False
    hook = hookfile.read()
    if python and hook == precommitbx_template(python):
        return repo_precommit_installed
    if "generated by precommitbx" in hook:
        return repo_precommit_legacy
    if hook:
        return repo_precommit_conflict
    return False


def clean_precommit_home():
    if precommit_home.check():
        print("Cleaning out existing pre-commit installation")
        precommit_home.remove(rec=1)


def python_version(python):
    result = clean_run([python, "--version"])
    if result.returncode:
        return False
    return result.stdout.strip().split(" ")[1]


def precommit_version(python):
    result = clean_run(
        [python, "-c", "import _precommitbx;print(_precommitbx.__version__)"],
        working_directory=precommit_home.join("..").join(".."),
    )
    if result.returncode:
        return False
    return result.stdout.strip()


def download(name, url, size, target, step=None):
    with tqdm(
        bar_format="{l_bar}{bar}| [{remaining}, {rate_fmt}]",
        desc="{:<18}".format(name),
        leave=False,
        total=size,
        unit="B",
        unit_scale=True,
    ) as bar:
        with target.open("wb", ensure=True) as fh:
            url_request = Request(url)
            with contextlib.closing(urlopen(url_request)) as socket:
                while True:
                    block = socket.read(4096)
                    if not block:
                        break
                    fh.write(block)
                    bar.update(len(block))
                    if step:
                        ProgressOverall.done_part(step, bar.n / size)
        if bar.n != size:
            raise RuntimeError(
                "Error downloading %s: received %d bytes, expected %d"
                % (url, bar.n, size)
            )


def download_openssl():
    archive = precommit_home / "openssl-{}.tar.gz".format(openssl_source_version)
    if archive.check() and archive.size() == openssl_source_size:
        return archive
    url = "https://www.openssl.org/source/openssl-{}.tar.gz".format(
        openssl_source_version
    )
    download("Downloading OpenSSL", url, openssl_source_size, archive, step="OpenSSL")
    return archive


def install_openssl():
    markerfile = precommit_home.join(".valid.openssl")
    if markerfile.check():
        return
    sourcedir = precommit_home / "openssl-{}".format(openssl_source_version)
    targz = download_openssl()
    ProgressOverall.done_step("OpenSSL")
    if sys.platform == "darwin":
        build_environment = {"KERNEL_BITS": "64"}
    else:
        build_environment = {}

    with Progress("Unpacking OpenSSL", 2416, step="OpenSSL") as bar:
        clean_run(
            ["tar", "xvfz", targz],
            working_directory=precommit_home,
            callback_stdout=bar.increment,
            stop_on_error="Error unpacking OpenSSL sources",
        )
    with Progress("Configuring OpenSSL", 351, step="OpenSSL") as bar:
        clean_run(
            [
                sourcedir.join("config"),
                "--prefix=%s" % precommit_home,
                "-fPIC",
                "no-hw",
            ],
            callback_stdout=bar.increment,
            environment_override=build_environment,
            stop_on_error="Error configuring OpenSSL sources",
            working_directory=sourcedir,
        )
    with Progress("Building OpenSSL", 1265, step="OpenSSL") as bar:
        clean_run(
            ["make"],
            callback_stdout=bar.increment,
            environment_override=build_environment,
            stop_on_error="Error building OpenSSL",
            working_directory=sourcedir,
        )
    with Progress("Installing OpenSSL", 2152, step="OpenSSL") as bar:
        clean_run(
            ["make", "install"],
            callback_stdout=bar.increment,
            environment_override=build_environment,
            stop_on_error="Error installing OpenSSL",
            working_directory=sourcedir,
        )
    markerfile.ensure()


def download_libffi():
    archive = precommit_home / "libffi-{}.tar.gz".format(libffi_source_version)
    if archive.check() and archive.size() == libffi_source_size:
        return archive
    url = "ftp://sourceware.org/pub/libffi/libffi-{}.tar.gz".format(
        libffi_source_version
    )
    download("Downloading libffi", url, libffi_source_size, archive, step="libffi")
    return archive


def install_libffi():
    markerfile = precommit_home.join(".valid.libffi")
    if markerfile.check():
        return
    sourcedir = precommit_home / "libffi-{}".format(libffi_source_version)
    targz = download_libffi()
    ProgressOverall.done_step("libffi")
    with Progress("Unpacking libffi", 358, step="libffi") as bar:
        clean_run(
            ["tar", "xvfz", targz],
            working_directory=precommit_home,
            callback_stdout=bar.increment,
            stop_on_error="Error unpacking libffi sources",
        )
    with Progress("Configuring libffi", 156, step="libffi") as bar:
        clean_run(
            [sourcedir.join("configure"), "--prefix=%s" % precommit_home],
            callback_stdout=bar.increment,
            stop_on_error="Error configuring libffi",
            working_directory=sourcedir,
        )
    with Progress("Building libffi", 76, step="libffi") as bar:
        clean_run(
            ["make"],
            callback_stdout=bar.increment,
            stop_on_error="Error building libffi",
            working_directory=sourcedir,
        )
    with Progress("Installing libffi", 63, step="libffi") as bar:
        clean_run(
            ["make", "install"],
            callback_stdout=bar.increment,
            stop_on_error="Error installing libffi",
            working_directory=sourcedir,
        )
    # Selectively copy output files to locations where Python setup can pick them up
    ffi_include = (
        precommit_home / "lib" / "libffi-{}".format(libffi_source_version) / "include"
    )
    if ffi_include.check():
        for f in ffi_include.listdir():
            f.copy(precommit_home / "include")
    ffi_lib = precommit_home.join("lib64")
    if ffi_lib.check():
        for f in precommit_home.join("lib64").listdir():
            f.copy(precommit_home / "lib")
    markerfile.ensure()


def download_sqlite():
    archive = precommit_home / "sqlite-autoconf-{}.tgz".format(sqlite_source_version)
    if archive.check() and archive.size() == sqlite_source_size:
        return archive
    url = "https://www.sqlite.org/2019/sqlite-autoconf-{}.tar.gz".format(
        sqlite_source_version
    )
    download("Downloading SQLite", url, sqlite_source_size, archive, step="SQLite")
    return archive


def install_sqlite():
    markerfile = precommit_home.join(".valid.sqlite")
    if markerfile.check():
        return
    sourcedir = precommit_home / "sqlite-autoconf-{}".format(sqlite_source_version)
    targz = download_sqlite()
    ProgressOverall.done_step("SQLite")

    with Progress("Unpacking SQLite", 44, step="SQLite") as bar:
        clean_run(
            ["tar", "xvfz", targz],
            working_directory=precommit_home,
            callback_stdout=bar.increment,
            stop_on_error="Error unpacking SQLite sources",
        )
    with Progress("Configuring SQLite", 116, step="SQLite") as bar:
        clean_run(
            [sourcedir.join("configure"), "--prefix=%s" % precommit_home],
            callback_stdout=bar.increment,
            stop_on_error="Error configuring SQLite sources",
            working_directory=sourcedir,
        )
    with Progress("Building SQLite", 17, step="SQLite") as bar:
        clean_run(
            ["make"],
            callback_stdout=bar.increment,
            stop_on_error="Error building SQLite",
            working_directory=sourcedir,
        )
    with Progress("Installing SQLite", 39, step="SQLite") as bar:
        clean_run(
            ["make", "install"],
            callback_stdout=bar.increment,
            stop_on_error="Error installing SQLite",
            working_directory=sourcedir,
        )
    markerfile.ensure()


def download_python():
    archive = precommit_home / "Python-{}.tgz".format(python_source_version)
    if archive.check() and archive.size() == python_source_size:
        return archive
    url = "https://www.python.org/ftp/python/{0}/Python-{0}.tgz".format(
        python_source_version
    )
    download("Downloading Python", url, python_source_size, archive, step="Python")
    return archive


def unpack_python():
    targz = download_python()
    ProgressOverall.done_step("Python")
    with Progress("Unpacking Python", 4183, step="Python") as bar:
        clean_run(
            ["tar", "xvfz", targz],
            working_directory=precommit_home,
            callback_stdout=bar.increment,
            stop_on_error="Error unpacking Python sources",
        )


def install_python(check_only=False):
    python3 = precommit_home.join("bin").join("python3")
    if python3.check():
        return python3
    if check_only:
        return False
    sourcedir = precommit_home / "Python-{}".format(python_source_version)
    ProgressOverall.start()
    if os.getenv("PARALLEL"):
        pool = multiprocessing.pool.ThreadPool(processes=4)
        pool.map(
            (lambda fn: fn()),
            (install_openssl, install_libffi, install_sqlite, unpack_python),
        )
        pool.close()
        ProgressOverall.done_package("OpenSSL")
        ProgressOverall.done_package("libffi")
        ProgressOverall.done_package("SQLite")
    else:
        install_openssl()
        ProgressOverall.done_package("OpenSSL")
        install_libffi()
        ProgressOverall.done_package("libffi")
        install_sqlite()
        ProgressOverall.done_package("SQLite")
        unpack_python()

    build_environment = {
        "LD_RUN_PATH": precommit_home.join("lib"),
        "LDFLAGS": "-L" + precommit_home.join("lib").strpath,
        "CPPFLAGS": "-I" + precommit_home.join("include").strpath,
    }

    with Progress("Configuring Python", 716, step="Python") as bar:
        clean_run(
            [
                sourcedir.join("configure"),
                "--prefix=%s" % precommit_home,
                "--with-openssl=%s" % precommit_home,
            ],
            callback_stdout=bar.increment,
            environment_override=build_environment,
            stop_on_error="Error configuring Python",
            working_directory=sourcedir,
        )
    try:
        parallel = str(multiprocessing.cpu_count())
    except NotImplementedError:
        parallel = "2"
    with Progress("Building Python", 461, step="Python") as bar:
        clean_run(
            ["make", "-j", parallel],
            callback_stdout=bar.increment,
            environment_override=build_environment,
            stop_on_error="Error building Python",
            working_directory=sourcedir,
        )
    compiled_python = sourcedir.join("python")
    if sys.platform == "darwin":
        compiled_python = sourcedir.join("python.exe")
    if not compiled_python.check():
        # in a parallel build 'make' might terminate before the build is complete
        for _ in trange(
            100,
            desc="Waiting for build results to appear...",
            bar_format="{l_bar}{bar}| [{remaining}]",
        ):
            if compiled_python.check():
                break
            time.sleep(0.1)
        else:
            sys.exit("Did not find build results")
    result = clean_run(
        [compiled_python, "-c", "import ssl; print(ssl._OPENSSL_API_VERSION)"],
        working_directory=sourcedir,
        stop_on_error="Python is missing SSL support",
    )
    if result.stdout.strip() != openssl_API_version:
        sys.exit("Python has not picked up correct OpenSSL headers")
    result = clean_run(
        [compiled_python, "-c", "import sqlite3; print(sqlite3.sqlite_version)"],
        working_directory=sourcedir,
        stop_on_error="Python is missing SQLite support",
    )
    if result.stdout.strip() != sqlite_version:
        sys.exit("Python has not picked up correct SQLite headers")
    result = clean_run(
        [compiled_python, "-c", "import ctypes"],
        working_directory=sourcedir,
        stop_on_error="Python is missing FFI support",
    )
    with Progress("Installing Python", 7770, step="Python") as bar:
        clean_run(
            ["make", "install"],
            working_directory=sourcedir,
            callback_stdout=bar.increment,
            stop_on_error="Error installing Python",
        )
    return python3


def install_precommit(python):
    ProgressOverall.start()
    for p in ("OpenSSL", "libffi", "SQLite", "Python"):
        ProgressOverall.done_package(p)
    with Progress("Installing Precommitbx", 28, step="Precommitbx") as bar:
        clean_run(
            [python, "-m", "pip", "install", py.path.local(__file__).dirname],
            callback_stdout=bar.increment,
            stop_on_error="Error installing precommitbx",
            working_directory=precommit_home,
        )


def list_all_repository_candidates():
    repositories = {}
    for module in sorted(libtbx.env.module_dict):
        module_paths = [
            py.path.local(abs(path))
            for path in libtbx.env.module_dict[module].dist_paths
            if path and (path / ".git").exists()
        ]
        if not module_paths:
            continue
        if len(module_paths) == 1:
            repositories[module] = module_paths[0]
        else:
            for path in module_paths:
                repositories[module + ":" + str(path)] = path
    import pkg_resources

    for ep in pkg_resources.iter_entry_points("libtbx.precommit"):
        path = py.path.local(ep.load().__path__[0])
        if path.join(".git").check():
            repositories[ep.name] = path
        elif path.dirpath().join(".git").check():
            repositories[ep.name] = path.dirpath()
    return repositories


def main():
    changes_required = False
    python = install_python(check_only=True)
    fix_things = "install" in sys.argv
    if python:
        py_ver = python_version(python)
    else:
        py_ver = False
    if py_ver != python_source_version and fix_things:
        python = install_python()
        py_ver = python_version(python)
    if py_ver == python_source_version:
        py_ver = GREEN + py_ver + NC
    elif py_ver:
        py_ver = YELLOW + py_ver + NC + " (expected: " + python_source_version + ")"
        changes_required = True
    else:
        py_ver = RED + "not installed" + NC
        changes_required = True
    if python:
        pc_ver = precommit_version(python)
        if pc_ver != precommitbx_version and fix_things:
            install_precommit(python)
            pc_ver = precommit_version(python)
        if pc_ver == precommitbx_version:
            pc_ver = GREEN + pc_ver + NC
        elif pc_ver:
            pc_ver = YELLOW + pc_ver + NC + " (expected: " + precommitbx_version + ")"
            changes_required = True
        else:
            pc_ver = RED + "not installed" + NC
            changes_required = True
        print("Precommit Python:", py_ver)
        print("Precommitbx:", pc_ver)

    print()
    print("Repositories:")
    repositories = list_all_repository_candidates()
    if "install" in sys.argv:
        paths = sys.argv[1:]
        paths.remove("install")
        for path in paths:
            path = py.path.local(".").join(path, abs=1)
            if path.basename in repositories:
                base = path.strpath
            else:
                base = path.basename
            repositories[base] = path
    for module in sorted(repositories):
        if not repositories[module].join(".pre-commit-config.yaml").check():
            print(repo_prefix.format(module), repo_no_precommit)
            continue
        message = (
            check_precommitbx_hook(repositories[module], python)
            or repo_precommit_available
        )
        if message != repo_precommit_installed and fix_things:
            install_precommitbx_hook(repositories[module], python)
            message = (
                check_precommitbx_hook(repositories[module], python)
                or repo_precommit_available
            )
        print(repo_prefix.format(module), message)
        if message != repo_precommit_installed:
            changes_required = True

    if changes_required:
        print()
        sys.exit(
            "To install pre-commit hooks run " + BOLD + "libtbx.precommit install" + NC
        )
