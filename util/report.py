#!/usr/bin/env python
#
# report.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function


class Array(object):
    """
    A class to represent an array

    """

    def __init__(self):
        """
        Initialise the array

        """
        self.name = ""
        self.title = ""
        self.data = None

    def as_dict(self):
        """
        Return as a dictionary

        :return: The dictionary

        """
        from collections import OrderedDict

        result = OrderedDict()
        result["title"] = self.title
        result["shape"] = self.data.all()
        result["data"] = list(self.data)
        return result

    def as_str(self, prefix=""):
        """
        Return as a string

        :return: The string

        """
        return ""


class Table(object):
    """
    A class to represent a table

    """

    def __init__(self):
        """
        Initialize the table

        """
        self.name = ""
        self.title = ""
        self.cols = []
        self.rows = []

    def as_dict(self):
        """
        Return as a dictionary

        :return: The dictionary

        """
        from collections import OrderedDict

        cols = OrderedDict()
        for col in self.cols:
            assert len(col) == 2
            cols[col[0]] = col[1]
        rows = []
        for row in self.rows:
            row_out = OrderedDict()
            for j, r in enumerate(row):
                row_out[self.cols[j][0]] = r
            rows.append(row_out)

        # Create the output
        result = OrderedDict()
        result["title"] = self.title
        result["cols"] = cols
        result["rows"] = rows
        return result

    def as_str(self, prefix=""):
        """
        Return the table as a string

        :return: The string

        """
        from libtbx.table_utils import format as table

        rows = [[col[1] for col in self.cols]]
        for i, row in enumerate(self.rows):
            rows.append([str(x) for x in row])
        text = [
            prefix + self.title,
            table(rows, has_header=True, justify="left", prefix=prefix),
            "",
        ]
        return "\n".join(text)


class Report(object):
    """
    A class to represent the report

    """

    def __init__(self):
        """
        Initialize the tables

        """
        self.tables = []
        self.arrays = []

    def add_array(self, array):
        """
        Add an array to the report

        :param array: The array to add

        """
        self.arrays.append(array)

    def add_table(self, table):
        """
        Add a table to the report

        :param table: The table to add

        """
        self.tables.append(table)

    def combine(self, other):
        """
        Combine two reports

        :param other: The other report

        """
        self.tables.extend(other.tables)
        self.arrays.extend(other.arrays)

    def as_dict(self):
        """
        Return the report as a dictionary

        :return: The dictionary

        """
        from collections import OrderedDict

        result = OrderedDict()
        result["tables"] = {table.name: table.as_dict() for table in self.tables}
        result["arrays"] = {array.name: array.as_dict() for array in self.arrays}
        return result

    def as_str(self, prefix=""):
        """
        Return the report as a string

        :return: The string

        """
        return "\n".join([table.as_str(prefix) for table in self.tables])

    def as_json(self):
        """
        Save the report as a json file

        :return: The json string

        """
        import json

        return json.dumps(self.as_dict(), indent=2)

    def as_xml(self):
        """
        Save the report as an xml file

        :return: The XML string

        """
        from xml.dom import minidom

        # Get the XML implementation
        impl = minidom.getDOMImplementation()

        # Create the XML document
        doc = impl.createDocument(None, None, None)

        # Create the document root
        root = doc.createElement("DIALS")

        # Function to process objects
        def process(root, obj):
            if isinstance(obj, dict):
                for key, value in obj.items():
                    root.appendChild(process(doc.createElement(key), value))
            elif isinstance(obj, list) or isinstance(obj, tuple):
                for i, value in enumerate(obj):
                    root.appendChild(process(doc.createElement("%d" % i), value))
            else:
                root.appendChild(doc.createTextNode(str(obj)))
            return root

        # Process the dictionary and append
        doc.appendChild(process(root, self.as_dict()))

        # Return the XML document
        return doc.toprettyxml(indent="  ")

    def as_file(self, filename):
        """
        Export as a file (either json or xml depending on extension

        :param filename: The filename

        """
        from os.path import splitext

        ext = splitext(filename)[1]
        with open(filename, "w") as outfile:
            if ext == ".json":
                outfile.write(self.as_json())
            elif ext == ".xml":
                outfile.write(self.as_xml())
            else:
                raise RuntimeError("Filename must be *.xml or *.json")
