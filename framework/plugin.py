#!/usr/bin/env python
#
# plugin.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from abc import ABCMeta, abstractmethod

def generate_phil_string(interface, extensions):

    choice_template = '''
    {name}
      .help = "{help}"
    {{
      algorithm = *{choices}
        .type = choice(multi={multi})
        .help = "Select the algorithm to use."
      {parameters}
    }}
    '''

    param_template = '''
      {name}
        .help = "{help}"
      {{
        {parameters}
      }}
    '''

    # Create all the phil parameters
    algorithms = []
    parameters = []
    for extension in extensions:
        algorithms.append(extension.name)
        try:
            parameters.append(param_template.format(
               name = extension.name,
               help = extension.__doc__,
               parameters = extension.phil))
        except Exception:
            pass

    # Get if multiple choice
    try:
        multi = interface.multi_choice
    except Exception:
        multi = False

    # Generate the choice string
    text = choice_template.format(
        name = interface.name,
        help = interface.__doc__,
        choices = ' '.join(algorithms),
        multi = multi,
        parameters = ' '.join(parameters))

    return text

def generate_single_phil(interface, extensions):
    from libtbx import phil
    return phil.parse(generate_phil_string(interface, extensions))


def generate_phil(interfaces):
    from libtbx import phil
    text = '\n'.join(generate_phil_string(*x) for x in interfaces.iteritems())
    return phil.parse(text)




class Factory(object):

    def __init__(self, plugins):
        self._plugins = plugins

    def create(self, name, *args, **kwargs):
        return self._plugins[name](*args, **kwargs)


def singleton(cls):
    instance = cls()
    instance.__call__ = lambda: instance
    return instance


@singleton
class Registry:

    def __init__(self):
        self._interfaces = set()

    def add(self, iface):
        self._interfaces.add(iface)

    def clear(self):
        self._interfaces.clear()

    def remove(self, iface):
        self._interfaces.remove(iface)

    def __len__(self):
        return len(self._interfaces)

    def __iter__(self):
        return iter(self._interfaces)

    def __contains__(self, iface):
        return iface in self._interfaces

    def extensions(self, cls):
        if cls not in self:
            raise TypeError('interface %s is not registered' % cls)

        stack = list(cls.__subclasses__())
        while len(stack) > 0:
            cls = stack.pop()
            yield cls
            stack.extend(cls.__subclasses__())

    def all_extensions(self):
        return dict((iface, list(self.extensions(iface))) for iface in self)

    def factory(self, iface):
        return Factory(dict((sc.name, sc) for sc in self.extensions(iface)))

    def interface_phil(self, iface):
        return generate_single_phil(iface, list(self.extensions(iface)))

    def global_phil(self):
        return generate_phil(self.all_extensions())


class Interface(ABCMeta):
    def __init__(self, name, bases, attrs):
        super(Interface, self).__init__(name, bases, attrs)
        if not hasattr(self, '_registered'):
            self._registered = True
            Registry.add(self)
