"""
Define an interface for observers and subjects.

A singleton decorator is also defined.
"""
from __future__ import absolute_import, division, print_function

import functools
from collections import OrderedDict


def singleton(cls):
    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]

    return getinstance


class Observer(object):
    def __init__(self):
        self.data = {}

    def update(self, subject):
        """Extract the relevant state information by inspecting the subject."""
        pass


class Subject(object):
    def __init__(self, events):
        self.observers = {}
        for event in events:
            self.observers[event] = OrderedDict()

    @staticmethod
    def notify_event(event):
        """Define a decorator for notifying of a specific event."""

        def wrap(method):
            @functools.wraps(method)
            def notify(self, *args, **kwargs):
                r = method(self, *args, **kwargs)
                self.notify(event)
                return r

            return notify

        return wrap

    def get_observers(self, event):
        return self.observers[event]

    def register_observer(self, event, observer, callback=None):
        if callback is None:
            callback = getattr(observer, "update")
        else:
            callback = getattr(observer, callback)
        self.get_observers(event)[observer] = callback

    def unregister_observer(self, event, observer):
        del self.get_observers(event)[observer]

    def notify(self, event):
        for callback in self.get_observers(event).values():
            callback(self)
