from __future__ import annotations

from dials.util import tabulate


class Array:
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

        return {
            "title": self.title,
            "shape": self.data.all(),
            "data": list(self.data),
        }

    def as_str(self, prefix=""):
        """
        Return as a string

        :return: The string
        """
        return ""


class Table:
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
        cols = {col[0]: col[1] for col in self.cols}
        rows = [{self.cols[j][0]: r for j, r in enumerate(row)} for row in self.rows]

        # Create the output
        return {
            "title": self.title,
            "cols": cols,
            "rows": rows,
        }

    def as_str(self, prefix=""):
        """
        Return the table as a string

        :return: The string
        """
        rows = [[col[1] for col in self.cols]]
        for i, row in enumerate(self.rows):
            rows.append([str(x) for x in row])
        text = [prefix + self.title, tabulate(rows, headers="firstrow"), ""]
        return "\n".join(text)


class Report:
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
        return {
            "tables": {table.name: table.as_dict() for table in self.tables},
            "arrays": {array.name: array.as_dict() for array in self.arrays},
        }

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
