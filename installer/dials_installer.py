from __future__ import absolute_import, division, print_function

import os
import shutil
import sys
import traceback

installer_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
libtbx_path = os.path.join(installer_path, "lib")
if libtbx_path not in sys.path:
    sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution


class installer(install_distribution.installer):
    organization = "dials"
    product_name = "DIALS"
    dest_dir_prefix = "dials"
    make_apps = []
    configure_modules = ["dials", "xia2", "iota", "prime"]
    include_gui_packages = True
    base_package_options = ["--dials"]
    installer_dir = installer_path
    modules = [
        # hot
        "annlib",
        "boost",
        "scons",
        "ccp4io",
        "msgpack-3.1.1",
        # base
        "cbflib",
        "cctbx_project",
        "gui_resources",
        "ccp4io_adaptbx",
        "annlib_adaptbx",
        "tntbx",
        "clipper",
        # dials
        "dxtbx",
        "dials",
        "xia2",
        "iota",
        "prime",
    ]
    flags = list(install_distribution.installer.flags)
    try:
        flags.remove("create_versioned_dispatchers")
    except ValueError:
        pass

    def product_specific_preinstallation_hook(self):
        prefix = os.path.abspath(self.options.prefix)
        if prefix.startswith(installer_path):
            sys.exit(
                "Invalid installation option: --prefix={givenprefix}\n\n"
                "Please install DIALS to a location outside of the installer directory.\n"
                "Suggested alternative: --prefix={suggestedprefix}".format(
                    givenprefix=self.options.prefix,
                    suggestedprefix=os.path.dirname(prefix),
                )
            )

    def reconfigure(self, log=None, *args, **kwargs):
        """Intercept any errors and print log excerpt"""
        try:
            return super(installer, self).reconfigure(log=log, *args, **kwargs)
        except Exception as e:
            if not self.options.verbose:
                print("\n" + " -=-" * 20)
                print("\nAn error occured during installation\n")
                print("Excerpt from installation log:")
                with open(log.name, "r") as fh:
                    for line in fh.readlines()[-30:]:
                        print(" :", line, end="")
                print("\nThis led to ", end="")
                sys.stdout.flush()
            traceback.print_exc()
            print("\n")
            sys.exit(
                "Please report this installation error to dials-support@lists.sourceforge.net"
            )

    def product_specific_prepackage_hook(self, directory):
        """
        Remove irrelevant files from installer.
        """
        self.print_header("Deflating installer")

        suffixes = ["B", "KB", "MB", "GB", "TB", "PB"]

        def humansize(nbytes):
            if nbytes == 0:
                return "0 B"
            i = 0
            while nbytes >= 1024 and i < len(suffixes) - 1:
                nbytes /= 1024.0
                i += 1
            f = ("%.2f" % nbytes).rstrip("0").rstrip(".")
            return "%s %s" % (f, suffixes[i])

        self._cleaned_size, self._cleaned_files = 0, 0

        def rmdir(subdir):
            fullpath = os.path.join(directory, subdir)
            if not os.path.exists(fullpath):
                print("Skipping", " " * 26, subdir)
                return
            num_files, total_size = 0, 0
            for dirpath, dirnames, filenames in os.walk(fullpath):
                for f in filenames:
                    fp = os.path.join(dirpath, f)
                    total_size += os.path.getsize(fp)
                    num_files += 1
            print(
                "Removing %9s, %4d files from %s"
                % (humansize(total_size), num_files, subdir)
            )
            shutil.rmtree(fullpath)
            self._cleaned_size = self._cleaned_size + total_size
            self._cleaned_files = self._cleaned_files + num_files

        def rmfile(filename):
            fullpath = os.path.join(directory, filename)
            if not os.path.exists(fullpath):
                print("Skipping", " " * 26, filename)
                return
            filesize = os.path.getsize(fullpath)
            print("Removing %9s, file %s" % (humansize(filesize), filename))
            os.remove(fullpath)
            self._cleaned_size = self._cleaned_size + filesize
            self._cleaned_files = self._cleaned_files + 1

        # Deduce matplotlib path
        # (base/lib/python2.??/site-packages/matplotlib-????/matplotlib)
        # (base/Python.framework/Versions/?.?/lib/python?.?/site-packages/matplotlib-(...) on MacOS)
        try:
            import matplotlib
            import inspect

            matplotpath = os.path.dirname(
                os.path.dirname(inspect.getsourcefile(matplotlib))
            )
            relpath = []
            matplotpath, d = os.path.split(matplotpath)
            relpath.append(d)
            while d and (d != "base"):
                matplotpath, d = os.path.split(matplotpath)
                relpath.append(d)
            if d == "base":
                relpath.reverse()
                # delete matplotlib tests
                matplotpath = os.path.join(*relpath)
                rmdir(os.path.join(matplotpath, "matplotlib", "tests"))
                rmdir(os.path.join(matplotpath, "mpl_toolkits", "tests"))

                # ...while we're here
                sitepath = os.path.dirname(matplotpath)
                rmdir(os.path.join(sitepath, "numpy/core/tests"))
                rmdir(os.path.join(sitepath, "numpy/doc"))
                rmdir(os.path.join(sitepath, "numpy/distutils/tests"))
                rmdir(os.path.join(sitepath, "numpy/f2py/docs"))

                pythonpath = os.path.dirname(sitepath)
                rmdir(os.path.join(pythonpath, "test"))
        except Exception:
            print("Could not deduce python package paths")

        rmdir("base/man")
        rmdir("base/share/doc")
        rmdir("base/share/gtk-doc")
        rmdir("base/share/hdf5_examples")
        rmdir("base/share/man")
        rmdir("build/dials_data")
        rmdir("build/precommitbx")
        rmdir("build/regression_data")
        rmdir("build/xia2_regression")
        for f in ("setpaths", "setpaths_debug", "setpaths_all", "unsetpaths"):
            for ext in (".sh", ".csh"):
                rmfile(os.path.join("build", f + ext))
        for p in (
            "chrono",
            "date_time",
            "detail",
            "filesystem",
            "program_options",
            "python",
            "system",
            "thread",
            "timer",
        ):
            rmdir(os.path.join("modules/boost/libs", p, "example"))
            rmdir(os.path.join("modules/boost/libs", p, "doc"))
            rmdir(os.path.join("modules/boost/libs", p, "test"))
            rmdir(os.path.join("modules/boost/libs", p, "tutorial"))
        rmdir("modules/boost/libs/date_time/xmldoc")
        rmdir("modules/cbflib/doc")
        rmdir("modules/cbflib/examples")
        rmdir("modules/cbflib/ply-3.2/doc")
        rmdir("modules/cbflib/ply-3.2/example")
        rmdir("modules/cbflib/ply-3.2/test")
        rmfile("modules/cbflib/idx-s00-20131106040304531.cbf")
        rmdir("modules/clipper/examples")
        print("-" * 60)
        print(
            "Deleted %d files, decrufting installation by %s\n"
            % (self._cleaned_files, humansize(self._cleaned_size))
        )


if __name__ == "__main__":
    installer(sys.argv[1:]).install()
