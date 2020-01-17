#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-

# Running bootstrap requires a minimum Python version of 2.7.

# To download this file:
# wget https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py
# or
# curl https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py > bootstrap.py

from __future__ import absolute_import, division, print_function

import argparse
import functools
import json
import multiprocessing
import os
import platform
import re
import shutil
import socket as pysocket
import stat
import subprocess
import sys
import tarfile
import time
import warnings
import zipfile

try:  # Python 3
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError
except ImportError:  # Python 2
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError

# Clean environment for subprocesses
clean_env = {
    key: value
    for key, value in os.environ.items()
    if key not in ("PYTHONPATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH")
}

# =============================================================================
# Locations for the files defining the conda environments
# Generally, they should reside in the repository of the main program
# defined by the builder. For example, the environment file for Phenix
# will be in the Phenix source tree.
conda_platform = {"Darwin": "osx-64", "Linux": "linux-64", "Windows": "win-64"}

devnull = open(os.devnull, "wb")  # to redirect unwanted subprocess output

allowed_ssh_connections = {}


def ssh_allowed_for_connection(connection):
    if connection not in allowed_ssh_connections:
        try:
            returncode = subprocess.call(
                ["ssh", "-oBatchMode=yes", "-T", connection],
                stdout=devnull,
                stderr=devnull,
            )
            # SSH errors lead to 255
            allowed_ssh_connections[connection] = returncode in (0, 1)
        except OSError:
            allowed_ssh_connections[connection] = False
    return allowed_ssh_connections[connection]


class conda_manager(object):
    def __init__(self):
        print()

        self.system = platform.system()

        # Find relevant conda base installation
        self.conda_base = os.path.realpath("miniconda")
        if os.name == "nt":
            self.conda_exe = os.path.join(self.conda_base, "Scripts", "conda.exe")
        else:
            self.conda_exe = os.path.join(self.conda_base, "bin", "conda")

        # default environment file for users
        self.environment_file = os.path.join(
            os.path.expanduser("~"), ".conda", "environments.txt"
        )

        if os.path.isdir(self.conda_base) and os.path.isfile(self.conda_exe):
            print("Using miniconda installation from", self.conda_base)
        else:
            print("Installing miniconda into", self.conda_base)
            self.install_miniconda(self.conda_base)

        # verify consistency and check conda version
        if not os.path.isfile(self.conda_exe):
            sys.exit("Conda executable not found at " + self.conda_exe)

        self.environments = self.update_environments()

        conda_info = json.loads(
            subprocess.check_output([self.conda_exe, "info", "--json"], env=clean_env)
        )
        if self.conda_base != os.path.realpath(conda_info["root_prefix"]):
            warnings.warn(
                "Expected conda base differs:",
                self.conda_base,
                "!=",
                os.path.realpath(conda_info["root_prefix"]),
            )
        for env in self.environments:
            if env not in conda_info["envs"]:
                print("Consistency check:", env, "not in environments:")
                print(conda_info["envs"])
                warnings.warn(
                    """
There is a mismatch between the conda settings in your home directory
and what "conda info" is reporting. This is not a fatal error, but if
an error is encountered, please check that your conda installation and
environments exist and are working.
""",
                    RuntimeWarning,
                )
        if conda_info["conda_version"] < "4.4":
            sys.exit(
                """
CCTBX programs require conda version 4.4 and greater to make use of the
common compilers provided by conda. Please update your version with
"conda update conda".
"""
            )

    def update_environments(self):
        """
    Read and check for existence of environment directories

    Returns
    -------
    environments: list
      List of paths that exist on the filesystem
    """
        environments = set()
        try:
            with open(self.environment_file) as f:
                paths = f.readlines()
            for env in paths:
                env = env.strip()
                if os.path.isdir(env):
                    environments.add(os.path.normpath(env))
        except IOError:
            pass

        env_dirs = [
            os.path.join(self.conda_base, "envs"),
            os.path.join(os.path.expanduser("~"), ".conda", "envs"),
        ]
        for env_dir in env_dirs:
            if os.path.isdir(env_dir):
                dirs = os.listdir(env_dir)
                for d in dirs:
                    d = os.path.join(env_dir, d)
                    if os.path.isdir(d):
                        environments.add(d)

        return environments

    def install_miniconda(self, location):
        """Download and install Miniconda3"""

        os_names = {"Darwin": "MacOSX", "Linux": "Linux", "Windows": "Windows"}
        filename = "Miniconda3-latest-{platform}-x86_64".format(
            platform=os_names[self.system]
        )
        if os.name == "nt":
            filename += ".exe"
        else:
            filename += ".sh"
        url_base = "https://repo.anaconda.com/miniconda/"
        url = url_base + filename
        filename = os.path.join(location, filename)

        print("Downloading {url}:".format(url=url), end=" ")
        result = download_to_file(url, filename)
        if result in (0, -1):
            sys.exit("Miniconda download failed")

        # run the installer
        if os.name == "nt":
            command = [
                filename,
                "/InstallationType=JustMe",
                "/RegisterPython=0",
                "/AddToPath=0",
                "/S",
                "/D=" + location,
            ]
        else:
            command = ["/bin/sh", filename, "-b", "-u", "-p", location]

        print()
        run_command(workdir=".", command=command, description="Installing Miniconda")

    def create_environment(self, python="36"):
        """
    Create the environment based on the builder and file. The
    environment name is "conda_base".

    Parameters
    ----------
    python: str
      If set, the specific Python version of the environment for the
      builder is used instead of the default. Current options are
      '27' and '36' for Python 2.7 and 3.6, respectively.
    """
        filename = os.path.join(
            "modules",
            "dials",
            ".conda-envs",
            "{builder}_py{version}_{platform}.txt".format(
                builder="dials",
                version=python,
                platform=conda_platform[platform.system()],
            ),
        )

        if not os.path.isfile(filename):
            raise RuntimeError(
                """The file {filename} is not available""".format(filename=filename)
            )

        # make a new environment directory
        prefix = os.path.realpath("conda_base")

        # install a new environment or update and existing one
        if prefix in self.environments:
            command = "install"
            text_messages = ["Updating", "update of"]
        else:
            command = "create"
            text_messages = ["Installing", "installation into"]
        command_list = [
            self.conda_exe,
            command,
            "--prefix",
            prefix,
            "--file",
            filename,
            "--yes",
            "--channel",
            "conda-forge",
            "--override-channels",
        ]
        if os.name == "nt":
            command_list = [
                "cmd.exe",
                "/C",
                " ".join(
                    [os.path.join(self.conda_base, "Scripts", "activate"), "base", "&&"]
                    + command_list
                ),
            ]
        print(
            "{text} {builder} environment with:\n  {filename}".format(
                text=text_messages[0], builder="dials", filename=filename
            )
        )
        for retry in range(5):
            retry += 1
            try:
                run_command(
                    workdir=".",
                    command=command_list,
                    description="Installing base directory",
                )
            except Exception:
                print(
                    """
*******************************************************************************
There was a failure in constructing the conda environment.
Attempt {retry} of 5 will start {retry} minute(s) from {t}.
*******************************************************************************
""".format(
                        retry=retry, t=time.asctime()
                    )
                )
                time.sleep(retry * 60)
            else:
                break
        else:
            sys.exit(
                """
The conda environment could not be constructed. Please check that there is a
working network connection for downloading conda packages.
"""
            )
        print(
            "Completed {text}:\n  {prefix}".format(text=text_messages[1], prefix=prefix)
        )
        with open(os.path.join(prefix, ".condarc"), "w") as fh:
            fh.write(
                """
channels:
  - conda-forge
"""
            )

        # on Windows, also download the Visual C++ 2008 Redistributable
        # use the same version as conda-forge
        # https://github.com/conda-forge/vs2008_runtime-feedstock
        if os.name == "nt":
            download_to_file(
                "https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x64.exe",
                os.path.join(prefix, "vcredist_x64.exe"),
            )

        # check that environment file is updated
        self.environments = self.update_environments()
        if prefix not in self.environments:
            raise RuntimeError(
                """
The newly installed environment cannot be found in
${HOME}/.conda/environments.txt.
"""
            )


_BUILD_DIR = "build"  # set by arg parser further on down


def tar_extract(workdir, archive):
    # using tarfile module rather than unix tar command which is not platform independent
    with tarfile.open(os.path.join(workdir, archive), errorlevel=2) as tar:
        tar.extractall(path=workdir)
        tarfoldername = os.path.join(
            workdir, os.path.commonprefix(tar.getnames()).split("/")[0]
        )
    # take full permissions on all extracted files
    module = os.path.join(workdir, tarfoldername)
    for root, dirs, files in os.walk(module):
        for fname in files:
            full_path = os.path.join(root, fname)
            os.chmod(
                full_path, stat.S_IREAD | stat.S_IWRITE | stat.S_IRGRP | stat.S_IROTH
            )


def run_command(command, workdir=_BUILD_DIR, description=None):
    print("===== Running in %s:" % workdir, description or " ".join(command))
    if workdir:
        try:
            os.makedirs(workdir)
        except OSError:
            pass
    try:
        p = subprocess.Popen(args=command, cwd=workdir, env=clean_env)
    except Exception as e:
        if isinstance(e, OSError):
            if e.errno == 2:
                executable = os.path.normpath(os.path.join(workdir, command[0]))
                raise RuntimeError("Could not run %s: File not found" % executable)
        if "child_traceback" in dir(e):
            print("Calling subprocess resulted in error; ", e.child_traceback)
        raise e

    try:
        p.wait()
    except KeyboardInterrupt:
        print("\nReceived CTRL+C, trying to terminate subprocess...\n")
        p.terminate()
        raise
    if p.returncode:
        sys.exit("Process failed with return code %s" % p.returncode)


def download_to_file(url, file, log=sys.stdout, status=True, cache=True):
    """Downloads a URL to file. Returns the file size.
       Returns -1 if the downloaded file size does not match the expected file
       size
       Returns -2 if the download is skipped due to the file at the URL not
       being newer than the local copy (identified by A. matching timestamp and
       size, or B. matching etag).
    """

    # Create directory structure if necessary
    if os.path.dirname(file):
        try:
            os.makedirs(os.path.dirname(file))
        except Exception:
            pass

    localcopy = os.path.isfile(file)

    # Get existing ETag, if present
    etag = None
    tagfile = "%s/.%s.etag" % os.path.split(os.path.abspath(file))
    if cache and os.path.isfile(tagfile):
        if not localcopy:
            # Having an ETag without a file is pointless
            os.remove(tagfile)
        else:
            tf = open(tagfile, "r")
            etag = tf.readline()
            tf.close()

    try:
        import ssl
        from ssl import SSLError
    except ImportError:
        ssl = None
        SSLError = None

    # Open connection to remote server
    try:
        if os.name == "nt" and "lbl.gov" in url:
            # Downloading from http://cci.lbl.gov/cctbx_dependencies caused
            # SSL: CERTIFICATE_VERIFY_FAILED error on Windows only as of today (why?).
            # Quick and dirty hack to disable ssl certificate verification.
            try:
                _create_unverified_https_context = ssl._create_unverified_context
            except AttributeError:
                # Legacy Python that doesn't verify HTTPS certificates by default
                pass
            except NameError:
                # ssl module was not loaded
                pass
            else:
                # Handle target environment that doesn't support HTTPS verification
                ssl._create_default_https_context = _create_unverified_https_context
        url_request = Request(url)
        if etag:
            url_request.add_header("If-None-Match", etag)
        if localcopy:
            # Shorten timeout to 7 seconds if a copy of the file is already present
            socket = urlopen(url_request, None, 7)
        else:
            socket = urlopen(url_request)
    except SSLError as e:
        # This could be a timeout
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            log.write("%s\n" % str(e))
            return -2
        # otherwise pass on the error message
        raise
    except (pysocket.timeout, HTTPError) as e:
        if isinstance(e, HTTPError) and etag and e.code == 304:
            # When using ETag. a 304 error means everything is fine
            log.write("local copy is current (etag)\n")
            return -2
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            log.write("%s\n" % str(e))
            return -2
        # otherwise pass on the error message
        raise
    except URLError as e:
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            log.write("%s\n" % str(e))
            return -2
        # if url fails to open, try using curl
        # temporary fix for old OpenSSL in system Python on macOS
        # https://github.com/cctbx/cctbx_project/issues/33
        command = ["/usr/bin/curl", "--http1.0", "-fLo", file, "--retry", "5", url]
        subprocess.call(command)
        socket = None  # prevent later socket code from being run
        try:
            received = os.path.getsize(file)
        except OSError:
            raise RuntimeError("Download failed")

    if socket is not None:
        try:
            file_size = int(socket.info().get("Content-Length"))
        except Exception:
            file_size = 0

        if os.path.isfile(tagfile):
            # ETag did not match, so delete any existing ETag.
            os.remove(tagfile)

        remote_mtime = 0
        try:
            remote_mtime = time.mktime(socket.info().getdate("last-modified"))
        except Exception:
            pass

        if file_size > 0:
            if remote_mtime > 0:
                # check if existing file matches remote size and timestamp
                try:
                    (
                        mode,
                        ino,
                        dev,
                        nlink,
                        uid,
                        gid,
                        size,
                        atime,
                        mtime,
                        ctime,
                    ) = os.stat(file)
                    if (size == file_size) and (remote_mtime == mtime):
                        log.write("local copy is current\n")
                        socket.close()
                        return -2
                except Exception:
                    # proceed with download if timestamp/size check fails for any reason
                    pass

            hr_size = (file_size, "B")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "kB")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "MB")
            log.write("%.1f %s\n" % hr_size)
            if status:
                log.write("    [0%")
                log.flush()

        received = 0
        block_size = 8192
        progress = 1
        # Allow for writing the file immediately so we can empty the buffer
        tmpfile = file + ".tmp"

        f = open(tmpfile, "wb")
        while True:
            block = socket.read(block_size)
            received += len(block)
            f.write(block)
            if status and (file_size > 0):
                while (100 * received / file_size) > progress:
                    progress += 1
                    if (progress % 20) == 0:
                        log.write("%d%%" % progress)
                    elif (progress % 2) == 0:
                        log.write(".")
                    log.flush()

            if not block:
                break
        f.close()
        socket.close()

        if status and (file_size > 0):
            log.write("]\n")
        else:
            log.write("%d kB\n" % (received / 1024))
        log.flush()

        # Do not overwrite file during the download. If a download temporarily fails we
        # may still have a clean, working (yet older) copy of the file.
        shutil.move(tmpfile, file)

        if (file_size > 0) and (file_size != received):
            return -1

        if remote_mtime > 0:
            # set file timestamp if timestamp information is available
            st = os.stat(file)
            atime = st[stat.ST_ATIME]  # current access time
            os.utime(file, (atime, remote_mtime))

        if cache and socket.info().get("ETag"):
            # If the server sent an ETAG, then keep it alongside the file
            open(tagfile, "w").write(socket.info().get("ETag"))

    return received


def unzip(archive, directory, trim_directory=0):
    """unzip a file into a directory."""
    print("===== Installing %s into %s" % (archive, directory))
    if not zipfile.is_zipfile(archive):
        raise Exception("%s is not a valid .zip file" % archive)
    z = zipfile.ZipFile(archive, "r")
    for member in z.infolist():
        is_directory = member.filename.endswith("/")
        filename = os.path.join(*member.filename.split("/")[trim_directory:])
        if filename != "":
            filename = os.path.normpath(filename)
            if "../" in filename:
                raise Exception(
                    "Archive %s contains invalid filename %s" % (archive, filename)
                )
            filename = os.path.join(directory, filename)
            upperdirs = os.path.dirname(filename)
            try:
                if is_directory and not os.path.exists(filename):
                    os.makedirs(filename)
                elif upperdirs and not os.path.exists(upperdirs):
                    os.makedirs(upperdirs)
            except Exception:
                pass
            if not is_directory:
                source = z.open(member)
                target = open(filename, "wb")
                shutil.copyfileobj(source, target)
                target.close()
                source.close()

                # Preserve executable permission, if set
                unix_executable = member.external_attr >> 16 & 0o111
                # rwxrwxrwx => --x--x--x => 0o111
                if unix_executable:
                    mode = os.stat(filename).st_mode
                    mode |= (mode & 0o444) >> 2  # copy R bits to X
                    # r--r--r-- => 0o444
                    os.chmod(filename, mode)
    z.close()


def set_git_repository_config_to_rebase(config):
    with open(config, "r") as fh:
        cfg = fh.readlines()

    branch, remote, rebase = False, False, False
    insertions = []
    for n, line in enumerate(cfg):
        if line.startswith("["):
            if branch and remote and not rebase:
                insertions.insert(0, (n, branch))
            if line.startswith("[branch"):
                branch = line.split('"')[1]
            else:
                branch = False
            remote, rebase = False, False
        if re.match(r"remote\s*=", line.strip()):
            remote = True
        if re.match(r"rebase\s*=", line.strip()):
            rebase = True
    if branch and remote and not rebase:
        insertions.insert(0, (n + 1, branch))
    for n, branch in insertions:
        print("  setting branch %s to rebase" % branch)
        cfg.insert(n, "\trebase = true\n")
    with open(config, "w") as fh:
        fh.write("".join(cfg))


def git(module, parameters, destination=None, reference=None):
    """Retrieve a git repository, either by running git directly
       or by downloading and unpacking an archive."""
    git_available = True
    try:
        subprocess.call(["git", "--version"], stdout=devnull, stderr=devnull)
    except OSError:
        git_available = False

    if destination is None:
        destination = os.path.join("modules", module)
    destpath, destdir = os.path.split(destination)

    if os.path.exists(destination):
        if git_available and os.path.exists(os.path.join(destination, ".git")):
            if (
                not open(os.path.join(destination, ".git", "HEAD"), "r")
                .read()
                .startswith("ref:")
            ):
                print(
                    "WARNING: Can not update existing git repository! You are not on a branch."
                )
                print(
                    "This may be legitimate when run eg. via Jenkins, but be aware that you cannot commit any changes"
                )
                return

            else:
                # This may fail for unclean trees and merge problems. In this case manual
                # user intervention will be required.
                # For the record, you can clean up the tree and *discard ALL changes* with
                #   git reset --hard origin/master
                #   git clean -dffx
                return run_command(
                    command=["git", "pull", "--rebase"], workdir=destination
                )

        print(
            "Existing non-git directory -- don't know what to do. skipping: %s" % module
        )
        return
    if isinstance(parameters, str):
        parameters = [parameters]
    git_parameters = []
    for source_candidate in parameters:
        if source_candidate.startswith("-"):
            git_parameters = source_candidate.split(" ")
            continue
        if not source_candidate.lower().startswith("http"):
            connection = source_candidate.split(":")[0]
            if not ssh_allowed_for_connection(connection):
                continue
        if source_candidate.lower().endswith(".git"):
            if not git_available:
                continue
            reference_parameters = []
            if reference is not None:
                if os.path.exists(reference) and os.path.exists(
                    os.path.join(reference, ".git")
                ):
                    reference_parameters = ["--reference", reference]
            cmd = (
                ["git", "clone", "--recursive"]
                + git_parameters
                + [source_candidate, destdir]
                + reference_parameters
                + ["--progress", "--verbose"]
            )
            run_command(command=cmd, workdir=destpath)
            if reference_parameters:
                # Sever the link between checked out and reference repository
                cmd = ["git", "repack", "-a", "-d"]
                run_command(command=cmd, workdir=destination)
                try:
                    os.remove(
                        os.path.join(
                            destination, ".git", "objects", "info", "alternates"
                        )
                    )
                except OSError:
                    set_git_repository_config_to_rebase(
                        os.path.join(destination, ".git", "config")
                    )
                    return 1
            set_git_repository_config_to_rebase(
                os.path.join(destination, ".git", "config")
            )
            # Show the hash for the checked out commit for debugging purposes
            run_command(command=["git", "rev-parse", "HEAD"], workdir=destination)
            return
        filename = "%s-%s" % (module, urlparse(source_candidate)[2].split("/")[-1])
        filename = os.path.join(destpath, filename)
        print("===== Downloading %s: " % source_candidate, end=" ")
        download_to_file(source_candidate, filename)
        unzip(filename, destination, trim_directory=1)
        return

    error = (
        "Cannot satisfy git dependency for module %s: None of the sources are available."
        % module
    )
    if not git_available:
        print(error)
        error = "A git installation has not been found."
    raise Exception(error)


def install_conda():
    conda_manager().create_environment()


def remove_files_by_extension(extension, workdir):
    cwd = os.getcwd()
    if workdir is not None:
        if os.path.exists(workdir):
            os.chdir(workdir)
        else:
            return
    print("\n  removing %s files in %s" % (extension, os.getcwd()))
    i = 0
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if name.endswith(extension):
                os.remove(os.path.join(root, name))
                i += 1
    os.chdir(cwd)
    print("  removed %d files" % i)


##### Modules #####
MODULES = {
    "scons": ["git", "-b 3.1.1", "https://github.com/SCons/scons/archive/3.1.1.zip"],
    "cctbx_project": [
        "git",
        "git@github.com:cctbx/cctbx_project.git",
        "https://github.com/cctbx/cctbx_project.git",
        "https://github.com/cctbx/cctbx_project/archive/master.zip",
    ],
    "boost": [
        "git",
        "git@github.com:cctbx/boost.git",
        "https://github.com/cctbx/boost.git",
        "https://github.com/cctbx/boost/archive/master.zip",
    ],
    "cbflib": [
        "git",
        "git@github.com:dials/cbflib.git",
        "https://github.com/dials/cbflib.git",
        "https://github.com/dials/cbflib/archive/master.zip",
    ],
    "annlib_adaptbx": [
        "git",
        "git@github.com:cctbx/annlib_adaptbx.git",
        "https://github.com/cctbx/annlib_adaptbx.git",
        "https://github.com/cctbx/annlib_adaptbx/archive/master.zip",
    ],
    "dials": [
        "git",
        "git@github.com:dials/dials.git",
        "https://github.com/dials/dials.git",
        "https://github.com/dials/dials/archive/master.zip",
    ],
    "dxtbx": [
        "git",
        "git@github.com:cctbx/dxtbx.git",
        "https://github.com/cctbx/dxtbx.git",
        "https://github.com/cctbx/dxtbx/archive/master.zip",
    ],
    "msgpack": [
        "curl",
        [
            "https://gitcdn.xyz/repo/dials/dependencies/dials-1.13/msgpack-3.1.1.tar.gz",
            "https://gitcdn.link/repo/dials/dependencies/dials-1.13/msgpack-3.1.1.tar.gz",
            "https://github.com/dials/dependencies/raw/dials-1.13/msgpack-3.1.1.tar.gz",
        ],
    ],
    "xia2": [
        "git",
        "git@github.com:xia2/xia2.git",
        "https://github.com/xia2/xia2.git",
        "https://github.com/xia2/xia2/archive/master.zip",
    ],
}
for module in (
    "annlib",
    "ccp4io",
    "ccp4io_adaptbx",
    "clipper",
    "eigen",
    "gui_resources",
    "tntbx",
):
    MODULES[module] = [
        "git",
        "git@github.com:dials/%s.git" % module,
        "https://github.com/dials/%s.git" % module,
        "https://github.com/dials/%s/archive/master.zip" % module,
    ]


###################################
##### Base Configuration      #####
###################################


class DIALSBuilder(object):
    # Configure these cctbx packages
    LIBTBX = [
        "cctbx",
        "cbflib",
        "dxtbx",
        "scitbx",
        "libtbx",
        "iotbx",
        "mmtbx",
        "smtbx",
        "gltbx",
        "wxtbx",
        "dials",
        "xia2",
        "prime",
        "iota",
        "--skip_phenix_dispatchers",
    ]

    def __init__(self, actions, options):
        """Create and add all the steps."""
        self.git_reference = options.git_reference
        self.steps = []
        # self.config_flags are only from the command line
        # LIBTBX can still be used to always set flags specific to a builder
        self.config_flags = options.config_flags or []

        # Add sources
        if "update" in actions:
            for m in sorted(MODULES):
                self.add_module(m)

        # always remove .pyc files
        self.remove_pyc()

        # Build base packages
        if "base" in actions:
            self.add_base()

        # Configure, make, get revision numbers
        if "build" in actions:
            self.add_configure()
            self.add_make()

        # Tests, tests
        if "tests" in actions:
            self.add_tests()

        if "build" in actions:
            self.add_refresh()

    def isPlatformMacOSX(self):
        return sys.platform.startswith("darwin")

    def remove_pyc(self):
        self.steps.append(
            functools.partial(remove_files_by_extension, ".pyc", "modules")
        )

    def run(self):
        for i in self.steps:
            i()

    def add_module(self, module):
        action = MODULES[module]
        method, parameters = action[0], action[1:]
        if len(parameters) == 1:
            parameters = parameters[0]
        if method == "curl":
            self._add_curl(module, parameters)
        elif method == "git":
            self._add_git(module, parameters)
        else:
            raise Exception("Unknown access method: %s %s" % (method, str(parameters)))

    def _add_download(self, url, to_file):
        if not isinstance(url, list):
            url = [url]

        def _download():
            for _url in url:
                for retry in (3, 3, 0):
                    print("===== Downloading %s: " % _url, end=" ")
                    try:
                        download_to_file(_url, to_file)
                        return
                    except Exception as e:
                        print("Download failed with", e)
                        if retry:
                            print("Retrying in %d seconds" % retry)
                            time.sleep(retry)
            raise RuntimeError("Could not download " + to_file)

        self.steps.append(_download)

    def _add_curl(self, module, url):
        if isinstance(url, list):
            filename = urlparse(url[0])[2].split("/")[-1]
        else:
            filename = urlparse(url)[2].split("/")[-1]
        # Google Drive URL does not contain the module name
        if filename == "uc":
            filename = module + ".gz"
        self._add_download(url, os.path.join("modules", filename))

        def _indirection():
            description = "extracting files from %s" % filename
            print("===== Running in modules:", description)
            tar_extract("modules", filename)

        self.steps.append(_indirection)

    def _add_git(self, module, parameters, destination=None):
        reference_repository_path = self.git_reference
        if reference_repository_path is None:
            if os.name == "posix" and pysocket.gethostname().endswith(".diamond.ac.uk"):
                reference_repository_path = (
                    "/dls/science/groups/scisoft/DIALS/repositories/git-reference"
                )
        if reference_repository_path:
            reference_repository_path = os.path.expanduser(
                os.path.join(reference_repository_path, module)
            )

        self.steps.append(
            functools.partial(
                git,
                module,
                parameters,
                destination=destination,
                reference=reference_repository_path,
            )
        )

    def _get_conda_python(self):
        """
    Helper function for determining the location of Python for the base
    and build actions.
    """
        if os.name == "nt":
            return os.path.join(os.getcwd(), "conda_base", "python.exe")
        elif self.isPlatformMacOSX():
            return os.path.join(
                "..", "conda_base", "python.app", "Contents", "MacOS", "python"
            )
        else:
            return os.path.join("..", "conda_base", "bin", "python")

    def add_command(self, command, description=None, workdir=None, args=None):
        if os.name == "nt":
            command = command + ".bat"
        # Relative path to workdir.
        workdir = workdir or [_BUILD_DIR]
        dots = [".."] * len(workdir)
        if workdir[0] == ".":
            dots = []
        if os.name == "nt":
            dots.extend([os.getcwd(), _BUILD_DIR, "bin", command])
        else:
            dots.extend([_BUILD_DIR, "bin", command])
        self.steps.append(
            functools.partial(
                run_command,
                command=[os.path.join(*dots)] + (args or []),
                description=description or command,
                workdir=os.path.join(*workdir),
            )
        )

    def add_refresh(self):
        self.add_command("libtbx.refresh", description="libtbx.refresh", workdir=["."])

    def add_base(self):
        self.steps.append(install_conda)

    def add_configure(self):
        if "--use_conda" not in self.config_flags:
            self.config_flags.append("--use_conda")

        configcmd = (
            [
                self._get_conda_python(),
                os.path.join(
                    "..", "modules", "cctbx_project", "libtbx", "configure.py"
                ),
            ]
            + self.LIBTBX
            + self.config_flags
        )
        self.steps.append(
            functools.partial(
                run_command,
                command=configcmd,
                description="run configure.py",
                workdir=_BUILD_DIR,
            )
        )

    def add_make(self):
        try:
            nproc = len(os.sched_getaffinity(0))
        except AttributeError:
            nproc = multiprocessing.cpu_count()
        self.add_command("libtbx.scons", args=["-j", str(nproc)])
        # run build again to make sure everything is built
        self.add_command("libtbx.scons", args=["-j", str(nproc)])

    def add_tests(self):
        self.add_command(
            "libtbx.pytest",
            args=["--regression", "-n", "auto"],
            description="test dxtbx",
            workdir=["modules", "dxtbx"],
        )
        self.add_command(
            "libtbx.pytest",
            args=["--regression", "-n", "auto"],
            description="test DIALS",
            workdir=["modules", "dials"],
        )


def run():
    prog = os.environ.get("LIBTBX_DISPATCHER_NAME")
    if prog is None or prog.startswith("python") or prog.endswith("python"):
        prog = os.path.basename(sys.argv[0])

    description = """
  You may specify one or more actions:
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Create conda environment
    build - Build
    tests - Run tests

  The default action is to run: update, base, build

  Example:

    python bootstrap.py update base build tests
  """

    parser = argparse.ArgumentParser(
        prog=prog,
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("action", nargs="*", help="Actions for building")
    parser.add_argument(
        "--git-reference",
        dest="git_reference",
        help="Path to a directory containing reference copies of repositories for faster checkouts.",
    )
    parser.add_argument(
        "--config-flags",
        dest="config_flags",
        help="""Pass flags to the configuration step. Flags should
be passed separately with quotes to avoid confusion (e.g
--config_flags="--build=debug" --config_flags="--enable_cxx11")""",
        action="append",
        default=[],
    )

    options = parser.parse_args()

    # Check actions
    allowed_actions = {"update", "base", "build", "tests"}
    actions = set(options.action or ["update", "base", "build"])
    unknown_actions = actions - allowed_actions
    if unknown_actions:
        sys.exit("Unknown action: %s" % ", ".join(unknown_actions))
    print("Performing actions:", " ".join(actions))
    DIALSBuilder(actions=actions, options=options).run()
    print("\nBootstrap success: %s" % ", ".join(actions))


if __name__ == "__main__":
    run()
