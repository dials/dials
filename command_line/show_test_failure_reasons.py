# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_test_failure_reasons

from __future__ import absolute_import, division, print_function

import re
import xml.etree.ElementTree as ET

try:
    tree = ET.parse("output.xml")
except IOError:
    exit(
        "generate output.xml with:\n"
        "pytest --regression -n auto --runxfail --junit-xml=output.xml"
    )
broken = tree.findall("testcase/failure") + tree.findall("testcase/error")
find_error_source = re.compile(
    r"\n((?:(?:/[^:\n]+/(?:modules|build)/[^:\n]+)|(?:[a-zA-Z][a-zA-Z0-9_./]*\.py)):[0-9]+):"
)

cause = {}
example = {}
for test in tree.iter("testcase"):
    for broken in test.findall("failure") + test.findall("error"):
        error_message = broken.attrib["message"]
        source = find_error_source.findall(broken.text)
        if source:
            error_message = source[-1]
        cause[error_message] = cause.setdefault(error_message, 0) + 1
        example[error_message] = {
            "source": "{t.attrib[file]}::{t.attrib[name]}".format(t=test),
            "text": broken.text,
        }

top_causes = sorted(((v, k) for k, v in cause.items()), reverse=True)


def filtered_output(output):
    lines = output.split("\n")
    error_lines = [n for n, line in enumerate(lines) if line.startswith(">")]
    if not error_lines:
        return ""
    return "\n".join(lines[error_lines[-1] :])


NC = "\033[0m"
RED = "\033[1;31m"
YELLOW = "\033[1;33m"

for cause in top_causes[:5]:
    print("{RED}{cause[0]}x {cause[1]}".format(RED=RED, cause=cause))
    print(
        "{YELLOW}pytest --regression --runxfail {source}{NC}".format(
            source=example[cause[1]]["source"], YELLOW=YELLOW, NC=NC
        )
    )
    print(filtered_output(example[cause[1]]["text"]))
    print()
