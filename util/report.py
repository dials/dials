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

from __future__ import division


class Table(object):
  '''
  A class to represent a table

  '''
  def __init__(self):
    '''
    Initialize the table

    '''
    self.title = ""
    self.cols = []
    self.rows = []

  def as_dict(self):
    '''
    Return as a dictionary

    :return: The dictionary

    '''
    from collections import OrderedDict
    rows = []
    for row in self.rows:
      row_out = OrderedDict()
      for j in range(len(row)):
        row_out[self.cols[j]] = row[j]
      rows.append(row_out)

    # Create the output
    result = OrderedDict()
    result['title'] = self.title
    result['cols'] = self.cols
    result['rows'] = rows
    return result


class Report(object):
  '''
  A class to represent the report

  '''

  def __init__(self):
    '''
    Initialize the tables

    '''
    self.tables = []

  def add_table(self, table):
    '''
    Add a table to the report

    :param table: The table to add

    '''
    self.tables.append(table)

  def combine(self, other):
    '''
    Combine two reports

    :param other: The other report

    '''
    self.tables.extend(other.tables)

  def as_dict(self):
    '''
    Return the report as a dictionary

    :return: The dictionary

    '''
    from collections import OrderedDict
    result = OrderedDict()
    result['tables'] = [ table.as_dict() for table in self.tables]
    return result

  def as_json(self):
    '''
    Save the report as a json file

    :return: The json string

    '''
    import json
    return json.dumps(self.as_dict(), indent=2)

  def as_xml(self):
    '''
    Save the report as an xml file

    :return: The XML string

    '''
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
        for key, value in obj.iteritems():
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
    '''
    Export as a file (either json or xml depending on extension

    :param filename: The filename

    '''
    from os.path import splitext
    ext = splitext(filename)[1]
    with open(filename, "w") as outfile:
      if ext == '.json':
        outfile.write(self.as_json())
      elif ext == '.xml':
        outfile.write(self.as_xml())
      else:
        raise RuntimeError('Filename must be *.xml or *.json')
