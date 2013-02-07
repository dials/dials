#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/call_method.hpp>
#include <iostream>
#include <memory>
#include <boost/noncopyable.hpp>

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

void callback_do_something(Base& a) {
    a.do_something();
}


template <typename BaseClass>
struct BaseClassInterface : BaseClass, boost::python::wrapper <BaseClass> {
    void do_something()
    {
        if (override f = this->get_override("do_something")) {
            f();
        } else {
            BaseClass::do_something();    
        }
    }

    void default_do_something() { 
        return this->BaseClass::do_something(); 
    }
};

struct BaseInterface : BaseClassInterface <Base> {};
struct DerivedInterface : BaseClassInterface <Derived> {};

template <typename T, typename B>
class_<T, boost::noncopyable> base_class_wrapper(const char *name) {
    return class_ <T, boost::noncopyable> (name)
        .def("do_something", 
            &B::do_something, 
            &T::default_do_something);
}

template <typename T, typename B>
class_<T, boost::noncopyable> derived_class_wrapper(const char *name) {
    return base_class_wrapper <T, B> (name);
}

void export_test_inheritance()
{
    class_<Base>("Base")
        .def("do_something", 
            &Base::do_something);
       
    class_<Derived, bases<Base> >("Derived");

    base_class_wrapper <BaseInterface, Base>("BaseInterface");
    derived_class_wrapper <DerivedInterface, Derived> ("DerivedInterface");


    def("callback_do_something", &callback_do_something);
};

}}}
