import os
import subprocess

# Flake8 validation
known_bad = {
    "algorithms/rs_mapper/__init__.py": {"F401", "F403"},
    "algorithms/shoebox/__init__.py": {"F401", "F403"},
    "nexus/__init__.py": {"F401", "F403"},
    "test/command_line/test_generate_distortion_maps.py": {"F841"},
}
failures = 0
try:
    flake8 = subprocess.run(
        [
            "flake8",
            "--exit-zero",
            "--max-line-length=88",
            "--select=E401,E711,E712,E713,E714,E721,E722,E901,F401,F402,F403,F405,F631,F632,F633,F811,F812,F821,F822,F841,F901,W191,W292,W293,W602,W603,W604,W605,W606",
        ],
        capture_output=True,
        check=True,
        encoding="latin-1",
        timeout=300,
    )
except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
    print(
        "##vso[task.logissue type=error;]flake8 validation failed with",
        str(e.__class__.__name__),
    )
    print(e.stdout)
    print(e.stderr)
    print("##vso[task.complete result=Failed;]flake8 validation failed")
    exit()
for line in flake8.stdout.split("\n"):
    if ":" not in line:
        continue
    filename, lineno, column, error = line.split(":", maxsplit=3)
    errcode, error = error.strip().split(" ", maxsplit=1)
    filename = os.path.normpath(filename)
    if errcode in known_bad.get(filename, {}):
        print("Ignoring warning", line)
    else:
        failures += 1
        print(
            f"##vso[task.logissue type=error;sourcepath={filename};"
            f"linenumber={lineno};columnnumber={column};code={errcode};]" + error
        )

if failures:
    print(f"##vso[task.logissue type=warning]Found {failures} flake8 violation(s)")
    print(f"##vso[task.complete result=Failed;]Found {failures} flake8 violation(s)")
