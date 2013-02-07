
from dials_util_ext import *

class PythonDerived (Base):

    def __init__(self):
        Base.__init__(self)

    def do_something(self):
        print "hello PythonDerived"

    def do_something_else(self):
        print "hello PythonDerived"

a = Base()
b = Derived()
c = ExtraDerived()
d = PythonDerived()

callback_do_something(a)
callback_do_something(b)
callback_do_something(c)
callback_do_something(d)
callback_do_something_else(c)

