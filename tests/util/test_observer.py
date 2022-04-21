"""
Test the observer module.
"""

from __future__ import annotations

from dials.util.observer import Observer, Subject, singleton


def test_singleton_decorator():
    """Test the singleton class decorator."""

    class TestClass:
        """Simplest class."""

        pass

    singleton_class = singleton(TestClass)

    c1 = singleton_class()
    c2 = singleton_class()

    assert c1 is c2

    c1 = TestClass()
    c2 = TestClass()

    assert c1 is not c2


def test_Observer():
    """Test for expected properties of the general observer."""

    observer = Observer()
    assert hasattr(observer, "data")
    observer.update([])


def test_Subject():
    """Test the subject class."""

    # First register events
    subject = Subject(events=["event1", "event2"])
    assert "event1" in subject.observers
    assert "event2" in subject.observers

    # Test the get_observers and register_observer and unregister_observer method
    observer = Observer()
    assert subject.get_observers("event1") == {}
    assert subject.get_observers("event2") == {}
    subject.register_observer("event1", observer)
    assert subject.get_observers("event1") == {observer: observer.update}
    assert subject.get_observers("event2") == {}
    subject.unregister_observer("event1", observer)
    assert subject.get_observers("event1") == {}
    assert subject.get_observers("event2") == {}

    # Try registering the observer with a different callback
    class MyObserver(Observer):
        """Test observer implementation class."""

        def __init__(self):
            self.count = 0

        def mycallback(self, _):
            """Custom callback method."""
            self.count += 1

    observer = MyObserver()
    subject = Subject(events=["event1", "event2"])
    subject.register_observer("event1", observer, callback="mycallback")
    assert subject.get_observers("event1") == {observer: observer.mycallback}
    assert subject.get_observers("event2") == {}

    # Now test call to notify.
    subject.notify("event1")
    assert observer.count == 1
    subject.notify("event2")
    assert observer.count == 1

    # Now test the notify decorator
    class MySubject(Subject):
        """Test subject implementation class."""

        def __init__(self):
            super().__init__(events=["event1"])

        @Subject.notify_event("event1")
        def test_method(self):
            """Test method for warpping."""
            return "test_return_value"

    mysubject = MySubject()
    myobserver = MyObserver()
    mysubject.register_observer("event1", myobserver, callback="mycallback")
    x = mysubject.test_method()
    assert myobserver.count == 1
    assert x == "test_return_value"
