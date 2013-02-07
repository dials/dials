#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/call_method.hpp>
#include <iostream>
#include <memory>

namespace dials { namespace util { namespace boost_python {

using namespace boost::python;

class Base {
public:
     virtual void do_something() { std::cout << "hello Base" << std::endl; }
};

class Derived : public Base {
public:
     virtual void do_something() { std::cout << "hello Derived" << std::endl; }
};

class ExtraDerived : public Derived {
public:
    virtual void do_something() { std::cout << "hello ExtraDerived" << std::endl; }
};

void callback_do_something(Base& a) {
    a.do_something();
}

struct BaseCallback : public Base {
public:

    BaseCallback(PyObject *p) : Base(), self(p) {}
    BaseCallback(PyObject *p, const Base& x) : Base(x), self(p) {}
    void do_something() { 
        call_method <void> (self, "do_something"); 
    }
    static void default_do_something(Base& self_) { 
        return self_.Base::do_something(); 
    }
 private:
    PyObject* self;
};

void export_test_inheritance()
{
    class_<Base, BaseCallback> ("Base")
        .def("do_something", &BaseCallback::default_do_something);
       
    class_<Derived, bases<Base> >("Derived")
        .def("do_something", &Derived::do_something);

    def("callback_do_something", &callback_do_something);
};

}}}
