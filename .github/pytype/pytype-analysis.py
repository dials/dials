from __future__ import annotations

import argparse
import collections
import itertools
import operator
import pathlib
import re

re_check = re.compile(r"^\[([0-9]+)/([0-9]+)\] check (.*)$")
re_import_error = re.compile(
    r"^File \"([^\"]+)\", line ([0-9]+), in [^:]*: Can't find module '([^']*)'. \[import-error\]$"
)
re_generic_error = re.compile(
    r"^File \"([^\"]+)\", line ([0-9]+), ([^\[]*)\[([^\]]+)\]$"
)
re_failed_line = re.compile(r"^FAILED: (.*)$")

parse_results = collections.namedtuple(
    "parse_results", "errors files missing_imports total_files"
)


def _github_formatting():
    log_format = collections.defaultdict(lambda: "%s")
    log_format["endgroup"] = "##[endgroup]%s"
    log_format["red"] = "\x1b[31;1m%s\x1b[0m"
    log_format["green"] = "\x1b[32;1m%s\x1b[0m"
    log_format["yellow"] = "\x1b[33;1m%s\x1b[0m"
    log_format["blue"] = "\x1b[34;1m%s\x1b[0m"
    log_format["cyan"] = "\x1b[36;1m%s\x1b[0m"
    log_format["group"] = "##[group]%s"
    return log_format


def _type_errors_were(number):
    if number == 1:
        return f"{number} type error was"
    else:
        return f"{number} type errors were"


def parse_options():
    parser = argparse.ArgumentParser(
        description="Parse pytype output for GitHub actions."
    )
    parser.add_argument(
        "pytypelog", metavar="LOG", type=str, help="pytype log to process"
    )
    parser.add_argument(
        "--base",
        dest="base",
        type=str,
        action="store",
        help="pytype log file to compare a run to",
    )
    parser.add_argument(
        "--github", action="store_true", help="Write log output in GitHub action format"
    )
    parser.add_argument(
        "--annotations",
        type=str,
        action="store",
        metavar="FILE",
        help="Write log for GitHub annotations",
    )
    return parser.parse_args()


def parse_pytypelog(logfile):
    logfile = pathlib.Path(logfile)
    if not logfile.is_file():
        exit("Could not find file " + str(logfile))

    errors = {}
    import_errors = collections.defaultdict(int)
    n_files = 0

    contents = [l.strip() for l in logfile.read_text().split("\n")]

    base_directory = None
    for line in contents:
        if re_failed_line.search(line):
            pyi_file_path = pathlib.Path(re_failed_line.search(line).group(1)).parts
            if ".pytype" in pyi_file_path:
                base_directory = pathlib.Path(
                    *pyi_file_path[: pyi_file_path.index(".pytype")]
                )
                break

    while contents:
        file_match = re_check.search(contents.pop(0))
        if not file_match:
            continue
        file_number, n_files, file_name = file_match.groups()

        # be optimistic
        errors[file_name] = []
        if not contents or re_check.search(contents[0]):
            continue
        while contents:
            # no more lines => end of entry
            if not contents[0]:
                # empty line => end of entry
                contents.pop(0)
                break
            if re_check.search(contents[0]):
                # next line = new match => end of entry
                break
            # any other output: something went wrong
            if re_import_error.search(contents[0]):
                import_errors[re_import_error.search(contents[0]).group(3)] += 1
            error_line = re_generic_error.search(contents[0])
            if error_line:
                if base_directory:
                    file_name_local = pathlib.Path(error_line.group(1)).relative_to(
                        base_directory
                    )
                    file_name_local = pathlib.Path(*file_name_local.parts[1:])
                else:
                    file_name_local = file_name  # best guess
                errors[file_name].append(
                    {
                        "filename_local": file_name_local,
                        "filename": error_line.group(1),
                        "line": error_line.group(2),
                        "description": error_line.group(3).strip(),
                        "error": error_line.group(4),
                    }
                )
            contents.pop(0)

    return parse_results(
        errors=sum(len(f) for f in errors.values()),
        missing_imports=import_errors,
        files=errors,
        total_files=n_files,
    )


if __name__ == "__main__":
    options = parse_options()

    if options.github:
        log_format = _github_formatting()
    else:
        log_format = collections.defaultdict(lambda: "%s")

    pytypelog = parse_pytypelog(options.pytypelog)
    print(f"{len(pytypelog.files)} out of {pytypelog.total_files} files were scanned")

    if pytypelog.missing_imports:
        print()
        print(
            log_format["group"]
            % (
                "Type information missing for %d packages:"
                % len(pytypelog.missing_imports)
            )
        )
        for package, count in sorted(pytypelog.missing_imports.items()):
            print(f"  {count:3d} reference(s) to {package}")
        print(log_format["endgroup"] % "")

    if pytypelog.errors:
        print()
        print(
            log_format["group"]
            % f"pytype found {pytypelog.errors} issues in the code base:"
        )
        error_types = collections.Counter(
            e["error"] for f in pytypelog.files.values() for e in f
        )
        for error, count in error_types.most_common():
            print(f"  {count:3d}x {error}")
        print(log_format["endgroup"] % "")
    else:
        print(log_format["green"] % "No errors found")

    if options.annotations and not options.base:
        with open(options.annotations, "w") as fh:
            for error in itertools.chain(*pytypelog.files.values()):
                fh.write(
                    "%s:%s:1:%s %s\n"
                    % (
                        str(error["filename_local"]),
                        error["line"],
                        error["error"],
                        error["description"],
                    )
                )

    if not options.base:
        exit(1 if pytypelog.errors else 0)

    base = parse_pytypelog(options.base)

    new_files = set(pytypelog.files) - set(base.files)
    if new_files:
        print()
        print(
            log_format["group"]
            % (f"{len(new_files)} modules were added or scanned for the first time:")
        )
        for f in sorted(new_files):
            print("  ", f)
        print(log_format["endgroup"] % "")

    def group_errors_by_type(errors):
        errors = sorted(errors, key=operator.itemgetter("error"))
        return {
            error_type: list(error_list)
            for error_type, error_list in itertools.groupby(
                errors, key=operator.itemgetter("error")
            )
        }

    annotations = []
    for f in set(pytypelog.files).intersection(base.files):
        past_errors = group_errors_by_type(base.files[f])
        current_errors = group_errors_by_type(pytypelog.files[f])
        differences = {
            error_type: len(current_errors.get(error_type, []))
            - len(past_errors.get(error_type, []))
            for error_type in (set(past_errors) | set(current_errors))
        }
        differences = {
            error_type: difference_count
            for error_type, difference_count in differences.items()
            if difference_count
        }
        if differences:
            fixed = sum(-d for d in differences.values() if d < 0)
            introduced = sum(d for d in differences.values() if d > 0)
            message = (
                (
                    log_format["green"] % (_type_errors_were(fixed) + " fixed ")
                    if fixed
                    else ""
                )
                + ("and " if fixed and introduced else "")
                + (
                    log_format["red"] % (_type_errors_were(introduced) + " introduced ")
                    if introduced
                    else ""
                )
                + f"in {f}"
            )

            print()
            print(message)
            for introduced_type in differences:
                if differences[introduced_type] > 0:
                    for error in current_errors[introduced_type]:
                        print(
                            f"   line {error['line']} {error['error']}: {error['description']}"
                        )
                        annotations.append(
                            "%s:%s:1:%s %s"
                            % (
                                str(error["filename_local"]),
                                error["line"],
                                error["error"],
                                error["description"],
                            )
                        )
    if options.annotations:
        with open(options.annotations, "w") as fh:
            fh.write("\n".join(annotations) + "\n")

    if annotations:
        exit(1)
