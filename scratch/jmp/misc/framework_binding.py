


class Method(object):

    def __call__(self, obj):
        print "hello"

    def __get__(self, obj, objtype):
        import types
        return types.MethodType(self, obj, objtype)

class A(object):
    pass


A.method = Method()

b = A()
b.method()
