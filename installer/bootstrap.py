#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-

# Running bootstrap requires a minimum Python version of 2.6.

# To download this file:
# wget https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py
# or
# curl https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py > bootstrap.py

from __future__ import absolute_import, division, print_function

import ntpath
import os
import os.path
import re
import shutil
import socket as pysocket
import stat
import subprocess
import sys
import tarfile
import time
import traceback

try:  # Python 3
    from urllib.parse import urlparse
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError
except ImportError:  # Python 2
    from urlparse import urlparse
    from urllib2 import urlopen, Request, HTTPError, URLError
import zipfile

try:
    import argparse
except ImportError:
    raise RuntimeError(
        """
The argparse module is required. If you are using Python 2.6, you may
need to install it separately. On CentOS 6, you can run

  yum install python-argpase

with administrative privileges.
"""
    )

_BUILD_DIR = "build"  # set by arg parser further on down

# Utility function to be executed on slave machine or called directly by standalone bootstrap script
def tar_extract(workdir, archive, modulename=None):
    try:
        # delete tar target folder if it exists
        if modulename and os.path.exists(modulename):

            def remShut(*args):
                func, path, _ = (
                    args
                )  # onerror returns a tuple containing function, path and exception info
                os.chmod(path, stat.S_IREAD | stat.S_IWRITE)
                os.remove(path)

            shutil.rmtree(modulename, onerror=remShut)
            # hack to work around possible race condition on Windows where deleted files may briefly
            # exist as phantoms and result in "access denied" error by subsequent IO operations
            cnt = 0
            while os.path.exists(modulename):
                time.sleep(1)
                cnt = cnt + 1
                if cnt > 5:
                    break
        # using tarfile module rather than unix tar command which is not platform independent
        tar = tarfile.open(os.path.join(workdir, archive), errorlevel=2)
        tar.extractall(path=workdir)
        tarfoldername = os.path.join(
            workdir, os.path.commonprefix(tar.getnames()).split("/")[0]
        )
        tar.close()
        # take full permissions on all extracted files
        module = os.path.join(workdir, tarfoldername)
        for root, dirs, files in os.walk(module):
            for fname in files:
                full_path = os.path.join(root, fname)
                os.chmod(
                    full_path,
                    stat.S_IREAD | stat.S_IWRITE | stat.S_IRGRP | stat.S_IROTH,
                )
        # rename to expected folder name, e.g. boost_hot -> boost
        # only rename if folder names differ
        if modulename:
            if modulename != tarfoldername:
                os.rename(tarfoldername, modulename)
    except Exception as e:
        raise Exception(
            "Extracting tar archive resulted in error: "
            + str(e)
            + "\n"
            + traceback.format_exc()
        )
        return 1
    return 0


# Mock commands to run standalone, without buildbot.
class ShellCommand(object):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def get_command(self):
        return self.kwargs["command"]

    def get_description(self):
        if "description" in self.kwargs:
            return self.kwargs["description"]
        return None

    def get_workdir(self):
        return self.kwargs.get("workdir", _BUILD_DIR)

    def get_environment(self):
        # gets environment from kwargs
        env = self.kwargs.get("env", None)
        if env:
            for key, item in env.items():
                if item is None:
                    env[key] = ""
                else:
                    env[key] = os.path.abspath(item)
            rc = os.environ
            rc.update(env)
            env = rc
        return env

    def run(self):
        command = self.get_command()
        description = self.get_description()
        workdir = self.get_workdir()
        env = self.get_environment()
        if not self.kwargs.get("quiet", False):
            if description:
                print("===== Running in %s:" % workdir, description)
            else:
                print("===== Running in %s:" % workdir, " ".join(command))
        if workdir:
            try:
                os.makedirs(workdir)
            except OSError:
                pass
        if command[0] == "tar":
            # don't think any builders explicitly calls tar but leave it here just in case
            modname = None
            if len(command) > 3 and command[3]:
                modname = command[3]
            return tar_extract(workdir, command[2], modname)
        if command[0] == "rm":
            # XXX use shutil rather than rm which is not platform independent
            for directory in command[2:]:
                if os.path.exists(directory):
                    print("Deleting directory : %s" % directory)
                    try:
                        shutil.rmtree(directory)
                    except OSError:
                        print("Strangely couldn't delete %s" % directory)
            return 0
        try:
            # if not os.path.isabs(command[0]):
            # executable path isn't located relative to workdir
            #  command[0] = os.path.join(workdir, command[0])
            stderr, stdout = None, None
            if self.kwargs.get("silent", False):
                stderr = stdout = open(os.devnull, "wb")
            p = subprocess.Popen(
                args=command, cwd=workdir, stdout=stdout, stderr=stderr, env=env
            )
        except Exception as e:  # error handling
            if not self.kwargs.get("haltOnFailure"):
                return 1
            if isinstance(e, OSError):
                if e.errno == 2:
                    executable = os.path.normpath(os.path.join(workdir, command[0]))
                    raise RuntimeError("Could not run %s: File not found" % executable)
            if "child_traceback" in dir(e):
                print("Calling subprocess resulted in error; ", e.child_traceback)
            raise e

        p.wait()
        if p.returncode != 0 and self.kwargs.get("haltOnFailure"):
            print("Process failed with return code %s" % (p.returncode))
            sys.exit(1)
        return p.returncode


class Toolbox(object):
    @staticmethod
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
            if sys.platform == "win32" and "lbl.gov" in url:
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
            subprocess.call(command, shell=False)
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
                from stat import ST_ATIME

                st = os.stat(file)
                atime = st[ST_ATIME]  # current access time
                os.utime(file, (atime, remote_mtime))

            if cache and socket.info().get("ETag"):
                # If the server sent an ETAG, then keep it alongside the file
                open(tagfile, "w").write(socket.info().get("ETag"))

        return received

    @staticmethod
    def unzip(archive, directory, trim_directory=0, verbose=False):
        """unzip a file into a directory."""
        if verbose:
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

    @staticmethod
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

    @staticmethod
    def git(
        module,
        parameters,
        destination=None,
        use_ssh=False,
        verbose=False,
        reference=None,
    ):
        """Retrieve a git repository, either by running git directly
       or by downloading and unpacking an archive."""
        git_available = True
        try:
            subprocess.call(
                ["git", "--version"],
                stdout=open(os.devnull, "wb"),
                stderr=open(os.devnull, "wb"),
            )
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
                    return ShellCommand(
                        command=["git", "pull", "--rebase"],
                        workdir=destination,
                        silent=False,
                        haltOnFailure=True,
                    ).run()

            print(
                "Existing non-git directory -- don't know what to do. skipping: %s"
                % module
            )
            return
        if isinstance(parameters, str):
            parameters = [parameters]
        git_parameters = []
        for source_candidate in parameters:
            if source_candidate.startswith("-"):
                git_parameters = source_candidate.split(" ")
                continue
            if not source_candidate.lower().startswith("http") and not use_ssh:
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
                )
                if verbose:
                    cmd = cmd + ["--progress", "--verbose"]
                returncode = ShellCommand(
                    command=cmd, workdir=destpath, silent=False
                ).run()
                if returncode:
                    return returncode  # no point trying to continue on error
                if reference_parameters:
                    # Sever the link between checked out and reference repository
                    cmd = ["git", "repack", "-a", "-d"]
                    returncode = ShellCommand(
                        command=cmd, workdir=destination, silent=False
                    ).run()
                    try:
                        os.remove(
                            os.path.join(
                                destination, ".git", "objects", "info", "alternates"
                            )
                        )
                    except OSError:
                        returncode = 1
                Toolbox.set_git_repository_config_to_rebase(
                    os.path.join(destination, ".git", "config")
                )
                if returncode:
                    return returncode  # no point trying to continue on error
                # Show the hash for the checked out commit for debugging purposes, ignore any failures.
                ShellCommand(
                    command=["git", "rev-parse", "HEAD"],
                    workdir=destination,
                    silent=False,
                ).run()
                return returncode
            filename = "%s-%s" % (module, urlparse(source_candidate)[2].split("/")[-1])
            filename = os.path.join(destpath, filename)
            if verbose:
                print("===== Downloading %s: " % source_candidate, end=" ")
            Toolbox.download_to_file(source_candidate, filename)
            Toolbox.unzip(filename, destination, trim_directory=1, verbose=verbose)
            return

        error = (
            "Cannot satisfy git dependency for module %s: None of the sources are available."
            % module
        )
        if not git_available:
            print(error)
            error = "A git installation has not been found."
        raise Exception(error)


class cleanup_ext_class(object):
    def __init__(self, filename_ext, workdir=None):
        self.filename_ext = filename_ext
        self.workdir = workdir

    def get_command(self):
        return "delete *%s in %s" % (self.filename_ext, self.workdir).split()

    def remove_ext_files(self):
        cwd = os.getcwd()
        if self.workdir is not None:
            if os.path.exists(self.workdir):
                os.chdir(self.workdir)
            else:
                return
        print("\n  removing %s files in %s" % (self.filename_ext, os.getcwd()))
        i = 0
        for root, dirs, files in os.walk(".", topdown=False):
            for name in files:
                if name.endswith(self.filename_ext):
                    os.remove(os.path.join(root, name))
                    i += 1
        os.chdir(cwd)
        print("  removed %d files" % i)

    def run(self):
        self.remove_ext_files()





##### Modules #####
class SourceModule(object):
    _modules = {}
    module = None
    anonymous = None

    def __init__(self):
        if not self._modules:
            self.update_subclasses()

    def items(self):
        return list(self._modules.items())

    @classmethod
    def update_subclasses(cls):
        for i in cls.__subclasses__():
            cls._modules[i.module] = i

    def get_module(self, module):
        if module in self._modules:
            return self._modules[module]
        raise KeyError("Unknown module: %s" % module)

    def get_url(self):
        repo = self.anonymous
        if not repo:
            raise Exception("No access method defined for module: %s" % self.module)
        return repo


# Core external repositories
# The trailing slashes ARE significant.
# These must all provide anonymous access.
class ccp4io_module(SourceModule):
    module = "ccp4io"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/ccp4io.gz",
            "https://drive.google.com/uc?id=1EF6AqowSrVnse7pRtRmIsvhS6Q0dsSLT&export=download",
        ],
    ]


class annlib_module(SourceModule):
    module = "annlib"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/annlib.gz",
            "https://drive.google.com/uc?id=1YD_KDXrfhJ5ryT97j4yxmbAPoecGLjg0&export=download",
        ],
    ]


class scons_module(SourceModule):
    module = "scons"
    anonymous = ["git", "-b 3.1.1", "https://github.com/SCons/scons/archive/3.1.1.zip"]


# Core CCTBX repositories
# These must all provide anonymous access.
class cctbx_module(SourceModule):
    module = "cctbx_project"
    anonymous = [
        "git",
        "git@github.com:cctbx/cctbx_project.git",
        "https://github.com/cctbx/cctbx_project.git",
        "https://github.com/cctbx/cctbx_project/archive/master.zip",
    ]


class boost_module(SourceModule):
    module = "boost"
    anonymous = [
        "git",
        "git@github.com:cctbx/boost.git",
        "https://github.com/cctbx/boost.git",
        "https://github.com/cctbx/boost/archive/master.zip",
    ]


class cbflib_module(SourceModule):
    module = "cbflib"
    anonymous = [
        "git",
        "git@github.com:dials/cbflib.git",
        "https://github.com/dials/cbflib.git",
        "https://github.com/dials/cbflib/archive/master.zip",
    ]


class ccp4io_adaptbx(SourceModule):
    module = "ccp4io_adaptbx"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/ccp4io_adaptbx.gz",
            "https://drive.google.com/uc?id=1X5kRE90KkV2yTEyF9zb-PHOjjRXjzYvx&export=download",
        ],
    ]


class annlib_adaptbx(SourceModule):
    module = "annlib_adaptbx"
    anonymous = [
        "git",
        "git@github.com:cctbx/annlib_adaptbx.git",
        "https://github.com/cctbx/annlib_adaptbx.git",
        "https://github.com/cctbx/annlib_adaptbx/archive/master.zip",
    ]


class tntbx_module(SourceModule):
    module = "tntbx"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/tntbx.gz",
            "https://drive.google.com/uc?id=1bDE_rF6iL0SeyplHSTNsfJyI1G1h7ZZv&export=download",
        ],
    ]


class clipper_module(SourceModule):
    module = "clipper"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/clipper.gz",
            "https://drive.google.com/uc?id=1xWAj59zoyVn26EoIuBrw7KLNRyGjS5wC&export=download",
        ],
    ]


class gui_resources_module(SourceModule):
    module = "gui_resources"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/gui_resources.gz",
            "https://drive.google.com/uc?id=1TTibOePamkUiIvwDJF-OMmdgX8jdgNUS&export=download",
        ],
    ]


class eigen_module(SourceModule):
    module = "eigen"
    anonymous = [
        "curl",
        [
            "http://cci.lbl.gov/repositories/eigen.gz",
            "https://drive.google.com/uc?id=138kErrF35WbnRRARqUczWaroao2w8p1A&export=download",
        ],
    ]


class dials_module(SourceModule):
    module = "dials"
    anonymous = [
        "git",
        "git@github.com:dials/dials.git",
        "https://github.com/dials/dials.git",
        "https://github.com/dials/dials/archive/master.zip",
    ]


class dxtbx_module(SourceModule):
    module = "dxtbx"
    anonymous = [
        "git",
        "git@github.com:cctbx/dxtbx.git",
        "https://github.com/cctbx/dxtbx.git",
        "https://github.com/cctbx/dxtbx/archive/master.zip",
    ]


class msgpack_module(SourceModule):
    module = "msgpack"
    anonymous = [
        "curl",
        [
            "https://gitcdn.xyz/repo/dials/dependencies/dials-1.13/msgpack-3.1.1.tar.gz",
            "https://gitcdn.link/repo/dials/dependencies/dials-1.13/msgpack-3.1.1.tar.gz",
            "https://github.com/dials/dependencies/raw/dials-1.13/msgpack-3.1.1.tar.gz",
        ],
    ]


class xia2_module(SourceModule):
    module = "xia2"
    anonymous = [
        "git",
        "git@github.com:xia2/xia2.git",
        "https://github.com/xia2/xia2.git",
        "https://github.com/xia2/xia2/archive/master.zip",
    ]


MODULES = SourceModule()

###################################
##### Base Configuration      #####
###################################


class DIALSBuilder(object):
    """Create buildbot configurations for CCI and CCTBX-like software."""

    # Checkout these codebases
    CODEBASES = [
        "boost",
        "cbflib",
        "cctbx_project",
        "dxtbx",
        "gui_resources",
        "ccp4io_adaptbx",
        "annlib_adaptbx",
        "tntbx",
        "clipper",
        "dials",
        "xia2",
    ]
    # Copy these sources from cci.lbl.gov
    HOT = ["annlib", "scons", "ccp4io", "eigen", "msgpack"]
    # Configure for these cctbx packages
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

    def __init__(
        self,
        python_base=None,
        hot=True,
        update=True,
        base=True,
        build=True,
        tests=True,
        auth=None,
        with_python=None,
        nproc=1,
        verbose=False,
        config_flags=[],
    ):
        if nproc is None:
            self.nproc = 1
        else:
            self.nproc = nproc
        """Create and add all the steps."""
        self.set_auth(auth)
        self.steps = []
        if self.isPlatformWindows():
            self.op = ntpath
        else:
            self.op = os.path
        self.name = "dials-dev"
        # Platform configuration.
        python_executable = "python3"
        if sys.platform == "win32":
            self.python_base = self.opjoin(
                *[os.getcwd(), "base", "bin", "python", python_executable]
            )
        else:
            self.python_base = self.opjoin(*["..", "base", "bin", python_executable])
        if with_python:
            self.python_base = with_python
        self.verbose = verbose
        # self.config_flags are only from the command line
        # LIBTBX can still be used to always set flags specific to a builder
        self.config_flags = config_flags

        # Add 'hot' sources
        if hot:
            list(map(self.add_module, self.HOT))

        # Add sources.
        if update:
            list(map(self.add_module, self.CODEBASES))

        # always remove .pyc files
        self.remove_pyc()

        # Build base packages
        if base:
            self.add_base()

        # Configure, make, get revision numbers
        if build:
            self.add_configure()
            self.add_make()
            self.add_install()

        # Tests, tests
        if tests:
            self.add_tests()

        if build:
            self.add_refresh()

    def isPlatformWindows(self):
        if sys.platform == "win32":
            return True
        return False

    def isPlatformMacOSX(self):
        if sys.platform.startswith("darwin"):
            return True
        return False

    def set_auth(self, auth):
        self.auth = auth or {}

    def remove_pyc(self):
        self.add_step(cleanup_ext_class(".pyc", "modules"))

    def shell(self, **kwargs):
        # Convenience for ShellCommand
        kwargs["haltOnFailure"] = kwargs.pop("haltOnFailure", True)
        kwargs["description"] = kwargs.get("description") or kwargs.get("name")
        kwargs["timeout"] = 60 * 60 * 2  # 2 hours
        if "workdir" in kwargs:
            kwargs["workdir"] = self.opjoin(*kwargs["workdir"])
        return ShellCommand(**kwargs)

    def run(self):
        for i in self.steps:
            i()

    def opjoin(self, *args):
        return self.op.join(*args)

    def add_step(self, step):
        """Add a step."""
        self.steps.append(step.run)

    def add_module(self, module, workdir=None, module_directory=None):
        action = MODULES.get_module(module)().get_url()
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

        class _download(object):
            def run(self):
                for _url in url:
                    for retry in (3, 3, 0):
                        print("===== Downloading %s: " % _url, end=" ")
                        try:
                            Toolbox().download_to_file(_url, to_file)
                            return
                        except Exception as e:
                            print("Download failed with", e)
                            if retry:
                                print("Retrying in %d seconds" % retry)
                                time.sleep(retry)
                raise RuntimeError("Could not download " + to_file)

        self.add_step(_download())

    def _add_curl(self, module, url):
        if isinstance(url, list):
            filename = urlparse(url[0])[2].split("/")[-1]
        else:
            filename = urlparse(url)[2].split("/")[-1]
        # Google Drive URL does not contain the module name
        if filename == "uc":
            filename = module + ".gz"
        self._add_download(url, os.path.join("modules", filename))
        self.add_step(
            self.shell(
                name="extracting files from %s" % filename,
                command=[
                    "python",
                    "-c",
                    "import sys; sys.path.append('..'); import bootstrap; \
       bootstrap.tar_extract('','%s')"
                    % filename,
                ],
                workdir=["modules"],
                description="extracting files from %s" % filename,
            )
        )

    def _add_git(self, module, parameters, destination=None):
        use_git_ssh = self.auth.get("git_ssh", False)
        reference_repository_path = self.auth.get("git_reference", None)
        if reference_repository_path is None:
            if os.name == "posix" and "diamond.ac.uk" in pysocket.gethostname():
                reference_repository_path = (
                    "/dls/science/groups/scisoft/DIALS/repositories/git-reference"
                )
        if reference_repository_path:
            reference_repository_path = os.path.expanduser(
                os.path.join(reference_repository_path, module)
            )

        class _indirection(object):
            def run(self):
                Toolbox().git(
                    module,
                    parameters,
                    destination=destination,
                    use_ssh=use_git_ssh,
                    verbose=True,
                    reference=reference_repository_path,
                )

        self.add_step(_indirection())

    def _get_conda_manager(self):
        """
    Helper function for determining the location of the conda environment
    """
        if __package__ is None:
            root_path = os.path.dirname(os.path.abspath(__file__))
            paths = [
                root_path,
                os.path.join(
                    root_path, "modules", "cctbx_project", "libtbx", "auto_build"
                ),
            ]
            for path in paths:
                if os.path.isfile(os.path.join(path, "install_conda.py")):
                    if path not in sys.path:
                        sys.path.append(path)
                        break
            from install_conda import conda_manager
        else:
            from .install_conda import conda_manager

        # drop output
        log = open(os.devnull, "w")

        # environment is provided, so do check that it exists
        check_file = False
        # base step has not run yet, so do not check if files exist
        conda_base = os.path.join("..", "conda_base")
        if self.isPlatformWindows():
            conda_base = os.path.join(os.getcwd(), "conda_base")
        # basic checks for python and conda
        m = conda_manager(
            root_dir=os.getcwd(), conda_env=conda_base, check_file=check_file, log=log
        )

        return m

    def _get_conda_python(self):
        """
    Helper function for determining the location of Python for the base
    and build actions.
    """
        try:
            m = self._get_conda_manager()
            return m.get_conda_python()
        except ImportError:  # modules directory is not available

            # -----------------------------------------------------------------------
            # duplicate logic from get_conda_python function in install_conda.py
            # since install_conda.py may not be available
            def m_get_conda_python(self):
                m_conda_python = os.path.join("bin", "python")
                if self.isPlatformWindows():
                    m_conda_python = self.op.join("python.exe")
                elif self.isPlatformMacOSX():
                    m_conda_python = os.path.join(
                        "python.app", "Contents", "MacOS", "python"
                    )
                return m_conda_python

            # -----------------------------------------------------------------------

            conda_python = None

            # (case 1)
            # use default location or file provided to --use-conda
            if True:
                conda_python = self.op.join(
                    "..", "conda_base", m_get_conda_python(self)
                )
                if self.isPlatformWindows():
                    conda_python = self.op.join(
                        os.getcwd(), "conda_base", m_get_conda_python(self)
                    )

            if conda_python is None:
                raise RuntimeError("A conda version of python could not be found.")

        return conda_python

    def add_command(self, command, name=None, workdir=None, args=None, **kwargs):
        if self.isPlatformWindows():
            command = command + ".bat"
        # Relative path to workdir.
        workdir = workdir or [_BUILD_DIR]
        dots = [".."] * len(workdir)
        if workdir[0] == ".":
            dots = []
        if sys.platform == "win32":
            dots.extend([os.getcwd(), _BUILD_DIR, "bin", command])
        else:
            dots.extend([_BUILD_DIR, "bin", command])
        self.add_step(
            self.shell(
                name=name or command,
                command=[self.opjoin(*dots)] + (args or []),
                workdir=workdir,
                **kwargs
            )
        )

    def add_test_command(
        self, command, name=None, workdir=None, args=None, haltOnFailure=False, **kwargs
    ):
        if name is None:
            name = "test %s" % command
        self.add_command(
            command,
            name=name,
            workdir=(workdir or ["tests", command]),
            args=args,
            haltOnFailure=haltOnFailure,
            **kwargs
        )

    def add_refresh(self):
        self.add_command("libtbx.refresh", name="libtbx.refresh", workdir=["."])

    def add_base(self):
        flags = ["--builder=dials"]
        # check for existing miniconda3 installation
        if not os.path.isdir("mc3"):
            flags.append("--install_conda")
        flags.append("--python=36")
        command = [
            "python",
            self.opjoin(
                "modules", "cctbx_project", "libtbx", "auto_build", "install_conda.py"
            ),
        ] + flags

        if len(command) > 0:
            print("Installing base packages using:\n  " + " ".join(command))
            self.add_step(self.shell(name="base", command=command, workdir=["."]))

    def add_configure(self):
        env = None

        if "--use_conda" not in self.config_flags:
            self.config_flags.append("--use_conda")
        self.python_base = self._get_conda_python()
        # conda python prefers no environment customizations
        # the get_environment function in ShellCommand updates the environment
        env = {"PYTHONPATH": None, "LD_LIBRARY_PATH": None, "DYLD_LIBRARY_PATH": None}

        configcmd = (
            [
                self.python_base,  # default to using our python rather than system python
                self.opjoin("..", "modules", "cctbx_project", "libtbx", "configure.py"),
            ]
            + self.LIBTBX
            + self.config_flags
        )
        self.add_step(
            self.shell(
                command=configcmd,
                workdir=[_BUILD_DIR],
                description="run configure.py",
                env=env,
            )
        )
        # Prepare saving configure.py command to file should user want to manually recompile Phenix
        fname = self.opjoin("config_modules.cmd")
        ldlibpath = ""
        confstr = ldlibpath + subprocess.list2cmdline(configcmd)
        if not self.isPlatformWindows():
            fname = self.opjoin("config_modules.sh")
            confstr = "#!/bin/sh\n\n" + confstr
        # klonky way of writing file later on, but it works
        self.add_step(
            self.shell(
                command=[
                    "python",
                    "-c",
                    'open(r"%s","w").write(r"""%s""" + "\\n")' % (fname, confstr),
                ],
                workdir=[_BUILD_DIR],
                description="save configure command",
            )
        )
        if not self.isPlatformWindows():
            self.add_step(
                self.shell(
                    command=["chmod", "+x", fname],
                    workdir=[_BUILD_DIR],
                    description="permit execution of config_modules.sh",
                )
            )

    def add_make(self):
        self.add_command("libtbx.scons", args=["-j", str(self.nproc)])
        # run build again to make sure everything is built
        self.add_command("libtbx.scons", args=["-j", str(self.nproc)])

    def add_install(self):
        """Run after compile, before tests."""
        self.add_command("mmtbx.rebuild_rotarama_cache", name="rebuild rotarama")

    def add_tests(self):
        self.add_test_command(
            "libtbx.pytest",
            args=["--regression", "-n", "auto"],
            workdir=["modules", "dxtbx"],
            haltOnFailure=True,
        )
        self.add_test_command(
            "libtbx.pytest",
            args=["--regression", "-n", "auto"],
            workdir=["modules", "dials"],
            haltOnFailure=True,
        )


def run():
    prog = os.environ.get("LIBTBX_DISPATCHER_NAME")
    if prog is None or prog.startswith("python") or prog.endswith("python"):
        prog = os.path.basename(sys.argv[0])

    description = """
  You may specify one or more actions:
    hot - Update static sources (scons, etc.)
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Build base dependencies (python, hdf5, wxWidgets, etc.)
    build - Build
    tests - Run tests

  The default action is to run: hot, update, base, build

  You can run the compilation step in parallel by providing a
  the number of processes using "--nproc".
  Complete build output is shown with "-v" or "--verbose".

  Finally, you may specify a specific Python interpreter
  using "--with-python".

  Example:

    python bootstrap.py hot update build tests
  """

    parser = argparse.ArgumentParser(
        prog=prog,
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("action", nargs="*", help="Actions for building")
    parser.add_argument(
        "--git-ssh",
        dest="git_ssh",
        action="store_true",
        help="Use ssh connections for git. This allows you to commit changes without changing remotes and use reference repositories.",
        default=False,
    )
    parser.add_argument(
        "--git-reference",
        dest="git_reference",
        help="Path to a directory containing reference copies of repositories for faster checkouts.",
    )
    parser.add_argument(
        "--with-python", dest="with_python", help="Use specified Python interpreter"
    )
    parser.add_argument("--nproc", help="number of parallel processes in compile step.")
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose output",
        default=False,
    )
    parser.add_argument(
        "--config-flags",
        "--config_flags",
        dest="config_flags",
        help="""Pass flags to the configuration step. Flags should
be passed separately with quotes to avoid confusion (e.g
--config_flags="--build=debug" --config_flags="--enable_cxx11")""",
        action="append",
        default=[],
    )

    parser.add_argument(
        "--build-dir",
        dest="build_dir",
        help="directory where the build will be. Should be at the same level as modules! default is 'build'",
        default="build",
        type=str,
    )

    options = parser.parse_args()
    args = options.action

    global _BUILD_DIR
    _BUILD_DIR = (
        options.build_dir
    )  # TODO: this is probably ok way to go with globalvar, but check and see

    # Check actions
    allowedargs = ["hot", "update", "base", "build", "tests"]
    args = args or ["hot", "update", "base", "build"]
    actions = []
    for arg in args:
        if arg not in allowedargs:
            raise ValueError("Unknown action: %s" % arg)
    for arg in allowedargs:
        if arg in args:
            actions.append(arg)

    print("Performing actions:", " ".join(actions))

    auth = {"git_ssh": options.git_ssh, "git_reference": options.git_reference}

    # Build
    DIALSBuilder(
        with_python=options.with_python,
        auth=auth,
        hot=("hot" in actions),
        update=("update" in actions),
        base=("base" in actions),
        build=("build" in actions),
        tests=("tests" in actions),
        nproc=options.nproc,
        verbose=options.verbose,
        config_flags=options.config_flags,
    ).run()
    print("\nBootstrap success: %s" % ", ".join(actions))


if __name__ == "__main__":
    run()
