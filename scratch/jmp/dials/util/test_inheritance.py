from __future__ import division

from dials_util_ext import *

class PythonDerived (DerivedInterface):

    def __init__(self):
        DerivedInterface.__init__(self)

    def do_something(self):
        print "hello PythonDerived"

class PythonDerived2 (PythonDerived):

    def __init__(self):
        PythonDerived.__init__(self)

    def do_something(self):
        print "hello PythonDerived2"


a = Base()
b = Derived()
d = PythonDerived()
e = PythonDerived2()

print "Callback Base"
callback_do_something(a)
print "Callback Derived"
callback_do_something(b)
print "Callback PythonDerived"
callback_do_something(d)

print "Do something Base"
a.do_something()
print "Do something Derived"
b.do_something()

print "Do something PythonDerived"
d.do_something()

print "Do something PythonDerived2"
e.do_something()

print "Callback Python derived2"
callback_do_something(e)
