
from dials_util_ext import *

class PythonDerived (Base):

    def __init__(self):
        Base.__init__(self)

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


callback_do_something(a)
callback_do_something(b)
callback_do_something(d)

a.do_something()
b.do_something()
d.do_something()
e.do_something()
callback_do_something(e)