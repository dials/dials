# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_tutorial_text

from __future__ import absolute_import, division, print_function

import glob
import os
import shlex
import shutil
import sys

try:
    import dials_data.download
except ImportError:
    dials_data = None
import libtbx.load_env  # required for libtbx.env.find_in_repositories
from libtbx.test_utils import open_tmp_directory
import procrunner


class Job(object):
    """Represents a step command to execute"""

    def __call__(self):
        """Run the command this job represents.

        Standard output is saved onto self.result, and then optionally filtered
        by the mangle_result function.

        :returns: Dictionary with keys 'cmd' and 'result'
        """
        print(self.cmd)
        self.result = Job.run_process(self.cmd)
        self.mangle_result()
        return {"cmd": self.cmd, "result": self.result["stdout"]}

    def mangle_result(self):
        """ function that can be overridden to change the return values after execution """
        pass

    @staticmethod
    def run_process(command):
        """Runs a command, prints running info and results the result, if success"""
        os.environ["DIALS_NOBANNER"] = "1"
        result = procrunner.run(shlex.split(command))
        print("running command took {:.2f} seconds\n".format(result["runtime"]))
        assert not result.returncode, "Command execution failed"
        return result


class JobWriter(object):
    """Tool to save job command and result to files in a fixed destination"""

    def __init__(self, directory):
        # Ensure this output path exists
        if not os.path.isdir(directory):
            os.makedirs(directory)

        self.directory = directory

    def __call__(self, cmd_filename, result_filename, job):
        """
        Save command and output to named files.

        :param cmd_filename:    Where to save a copy of the run command
        :param result_filenam:  Where to save the result output
        :param job:             The result dictionary from calling a Job
        """
        with open(os.path.join(self.directory, cmd_filename), "w") as f:
            f.write(job["cmd"])
        with open(os.path.join(self.directory, result_filename), "w") as f:
            f.write(job["result"])


class Processing_Tutorial(object):
    """Command steps for generating the logs for the processing tutorial"""

    class dials_import(Job):
        def __init__(self):
            # find i04 bag training data
            if not dials_data:
                raise RuntimeError(
                    "You need to install the dials_data python package first.\n"
                    "Run libtbx.pip install dials_data"
                )

            df = dials_data.download.DataFetcher()
            dataset = df("thaumatin_i04").join("th_8_2_0*cbf").strpath

            self.cmd = "dials.import {}".format(dataset)

    class dials_find_spots(Job):
        cmd = "dials.find_spots imported.expt nproc=4"

    class dials_index(Job):
        cmd = "dials.index imported.expt strong.refl"

    class dials_refine_bravais_settings(Job):
        cmd = "dials.refine_bravais_settings indexed.expt indexed.refl"

    class dials_reindex(Job):
        cmd = "dials.reindex indexed.refl change_of_basis_op=a,b,c"

    class dials_refine(Job):
        cmd = "dials.refine bravais_setting_9.expt reindexed.refl scan_varying=false"

    class dials_sv_refine(Job):
        cmd = "dials.refine refined.expt refined.refl scan_varying=true"

    class dials_integrate(Job):
        cmd = "dials.integrate refined.expt refined.refl nproc=4"

    class dials_report(Job):
        cmd = "dials.report integrated.expt integrated.refl"

        def mangle_result(self):
            self.result["stdout"] = open("dials-report.html").read()

    class dials_export(Job):
        cmd = "dials.export integrated.refl integrated.expt"


def generate_processing_detail_text_thaumatin():
    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="generate_tutorial_text")
    os.chdir(tmp_dir)

    try:
        import_job = Processing_Tutorial.dials_import()()

        find_spots_job = Processing_Tutorial.dials_find_spots()()

        index_job = Processing_Tutorial.dials_index()()

        refine_bravais_settings_job = (
            Processing_Tutorial.dials_refine_bravais_settings()()
        )

        reindex_job = Processing_Tutorial.dials_reindex()()

        refine_job = Processing_Tutorial.dials_refine()()

        sv_refine_job = Processing_Tutorial.dials_sv_refine()()

        integrate_job = Processing_Tutorial.dials_integrate()()

        report_html_job = Processing_Tutorial.dials_report()()

        export_job = Processing_Tutorial.dials_export()()

        # if we got this far, assume it is okay to overwrite the logs
        dials_dir = libtbx.env.find_in_repositories("dials")
        result_dir = os.path.join(
            dials_dir, "doc", "sphinx", "documentation", "tutorials", "logs"
        )

        job_writer = JobWriter(result_dir)
        job_writer("dials.import.cmd", "dials.import.log", import_job)
        job_writer("dials.find_spots.cmd", "dials.find_spots.log", find_spots_job)
        job_writer("dials.index.cmd", "dials.index.log", index_job)
        job_writer(
            "dials.refine_bravais_settings.cmd",
            "dials.refine_bravais_settings.log",
            refine_bravais_settings_job,
        )
        job_writer("dials.reindex.cmd", "dials.reindex.log", reindex_job)
        job_writer("dials.refine.cmd", "dials.refine.log", refine_job)
        job_writer("dials.sv_refine.cmd", "dials.sv_refine.log", sv_refine_job)
        job_writer("dials.integrate.cmd", "dials.integrate.log", integrate_job)
        job_writer("dials.report.cmd", "dials-report.html", report_html_job)
        job_writer("dials.export.cmd", "dials.export.log", export_job)

        print("Updated result files written to {}".format(result_dir))

    finally:
        os.chdir(cwd)
        # clean up tmp dir
        shutil.rmtree(tmp_dir)


def generate_processing_detail_text_betalactamase():
    """Generate the text for Beta-lactamase versions of detail processing tutorial"""

    # Move to a temporary directory for processing
    cwd = os.path.abspath(os.curdir)
    tmp_dir = open_tmp_directory(suffix="generate_tutorial_text_ccp4")
    os.chdir(tmp_dir)

    # Find/validate the data input - until we've decided to integrate this
    # into the main release, have a DLS default or otherwise let it be
    # specified via an environment variable.
    DATA_PATH = (
        "/dls/i03/data/2017/mx19576-1/tutorial_data/summed/summed/C2sum_1*.cbf.gz"
    )
    # If not on dls systems, look for an environment variable
    if not glob.glob(DATA_PATH):
        # Firstly, look for the old path
        if "CCP4_TUTORIAL_DATA" in os.environ:
            DATA_PATH = os.environ.get("CCP4_TUTORIAL_DATA")
        # Otherwise, look for the new path that we tell the user about
        elif "BETALACTAMASE_TUTORIAL_DATA" in os.environ:
            base_path = os.environ.get("BETALACTAMASE_TUTORIAL_DATA")
            DATA_PATH = os.path.join(base_path, "C2sum_1*.cbf.gz")

    if not DATA_PATH or not any(glob.glob(DATA_PATH)):
        sys.exit(
            """Error:  Could not find Betalactamase data: skipping text generation.
        Please download C2sum_1 from https://zenodo.org/record/1014387 and extract,
        then set environment variable BETALACTAMASE_TUTORIAL_DATA to the folder with C2sum_1_*.cbf.gz"""
        )

    print("Using data: {}".format(DATA_PATH))
    # Work out where we are writing the output files to; in-source
    dials_dir = libtbx.env.find_in_repositories("dials")
    OUTPUT_DIR = os.path.join(
        dials_dir,
        "doc",
        "sphinx",
        "documentation",
        "tutorials",
        "logs_detail_betalactamase",
    )
    # Ensure this output path exists
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    print("Writing output logs to {}".format(OUTPUT_DIR))

    # Make an ordered list of named steps and associated commands
    commands = [
        ("dials.import", "dials.import {}".format(DATA_PATH)),
        ("dials.find_spots", "dials.find_spots imported.expt nproc=4"),
        ("dials.index", "dials.index imported.expt strong.refl"),
        (
            "dials.refine_bravais_settings",
            "dials.refine_bravais_settings indexed.expt indexed.refl",
        ),
        ("dials.reindex", "dials.reindex indexed.refl change_of_basis_op=a+b,-a+b,c"),
        (
            "dials.refine",
            "dials.refine bravais_setting_2.expt reindexed.refl scan_varying=false",
        ),
        ("dials.sv_refine", "dials.refine refined.expt refined.refl scan_varying=true"),
        ("dials.integrate", "dials.integrate refined.expt refined.refl nproc=4"),
        ("dials.report", "dials.report integrated.expt integrated.refl"),
        ("dials.export", "dials.export integrated.refl integrated.expt"),
    ]

    job_writer = JobWriter(OUTPUT_DIR)

    # Protect against errors so we can tidy up afterwards
    try:
        # Run each step, and write an output log
        for name, command in commands:
            result = Job.run_process(command)["stdout"]
            # Write a copy of the command, and the output log
            job_writer(
                "{}.cmd".format(name),
                "{}.log".format(name),
                {"cmd": command, "result": result},
            )

        # Because it's hard to robustly extract the *final* unindexed count
        # using sphinx's built-in literalinclude, we extract it specially here
        extract_last_indexed_spot_count(os.path.join(OUTPUT_DIR, "dials.index.log"))

        # Report step is special; we want the dials-report.html file instead
        shutil.copy("dials-report.html", OUTPUT_DIR)

        print("Updated result files written to {}".format(OUTPUT_DIR))

    finally:
        # Remove our intermediatary files
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)


def find_in_line(string, lines, start=0):
    """Find the next line index containing a given string"""
    for n, line in enumerate(lines[start:], start):
        if string in line:
            return n
    raise RuntimeError("Could not find line '{}' in lines".format(string))


def write_extract(filename, start, end, lines):
    """Write lines to a file, in the correct line position, with markers.

    This can be used to provide sphinx with an easily extractable literalinclude
    block, that preserves the correct line numbers from the original file.
    """
    assert start > 0
    out_lines = []
    for n, line in enumerate(lines):
        if n == start - 1:
            out_lines.append("[START_EXTRACT]")
        elif n >= start and n <= end:
            out_lines.append(line.strip())
        elif n == end + 1:
            out_lines.append("[END_EXTRACT]")
        else:
            out_lines.append("")

    with open(filename, "wt") as f:
        f.write("\n".join(out_lines))


def extract_last_indexed_spot_count(path):
    """Extract the last 'unindexed' entry from an index log."""
    with open(path) as f:
        lines = f.readlines()

    # Find the last entry for '# unindexed'
    next_ui = len(lines) - find_in_line("# unindexed", list(reversed(lines))) - 1

    # We now have the last entry for '# unindexed'. Find the end of the table
    end_ui = find_in_line("---", lines, next_ui + 2)
    # Range to extract = (next_ui, end_ui)
    dest = os.path.dirname(path)
    write_extract(
        os.path.join(dest, "dials.index.log.extract_unindexed"),
        next_ui - 1,
        end_ui,
        lines,
    )


if __name__ == "__main__":
    if not dials_data:
        sys.exit(
            "You need to install the dials_data python package first.\n"
            "Run libtbx.pip install dials_data"
        )

    if "-h" in sys.argv or "--help" in sys.argv:
        print("Usage: dev.dials.generate_tutorial_text [--beta | --thaum]")
        sys.exit(0)

    # As a quick development hack, add option for only the newer process
    if "--beta" not in sys.argv:
        print("Generating thaumatin tutorial")
        generate_processing_detail_text_thaumatin()
    if "--thaum" not in sys.argv:
        print("Generating betalactamase tutorial")
        generate_processing_detail_text_betalactamase()
