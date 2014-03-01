"""
Remote-Callable Method Decorator

This decorator is syntactic sugar for defining methods to be called
from the remote end of a rpc object.

See :class:`DecoratorTest` for an example of how the decorator works.

EXAMPLES::

    sage: from sage.rpc.core.decorator import DecoratorTest
    sage: test = DecoratorTest()
    Register foo to call _impl_foo
    sage: test
    <sage.rpc.core.decorator.DecoratorTest object at 0x...>

TESTS::

    sage: from sage.rpc.core.decorator import remote_callable, remote_callable_iter
    sage: class Test2(object):
    ....:     def __init__(self):
    ....:         list(remote_callable_iter(self))
    ....:     @remote_callable('foo')
    ....:     def _foo_impl_1(self):
    ....:         pass
    ....:     @remote_callable('foo')   # duplicate name is not allowed
    ....:     def _foo_impl_2(self):
    ....:         pass
    sage: Test2()
    Traceback (most recent call last):
    ...
    RuntimeError: RPC method name "foo" already in use
"""

import sys
import types


class remote_callable(object):
    
    def __init__(self, name):
        """
        Record the RPC method name for later use.
        
        INPUT:

        - ``name`` -- string. A valid name for a RPC method.

        EXAMPLES::

            sage: from sage.rpc.core.decorator import remote_callable
            sage: def f(x): return x
            sage: holder = remote_callable('call.f')(f)
            sage: holder.name
            'call.f'
            sage: holder.method
            <function f at 0x...>
        """
        self.name = name

    def __call__(self, method):
        """
        Use the decorator to record the RPC implementation for later use.
        
        INPUT:

        - ``method`` -- a bound method. The RPC implementation to call.

        EXAMPLES::

            sage: from sage.rpc.core.decorator import remote_callable
            sage: def f(x): return x
            sage: holder = remote_callable('call.f')(f)
            sage: holder.name
            'call.f'
            sage: holder.method
            <function f at 0x...>
        """
        self.method = method
        return self


def remote_callable_iter(instance):
    """
    Iterate over the remote callable methods

    The remote callable methods are the methods marked with the
    :class:`remote_callable` decorator. This function must be called
    on the instance during the initialization (i.e. from
    ``__init__``), where it replaces the decorated method with the
    actual method. While doing so, it yields the rpc method name
    defined in the decorator, so it can be registered with the rpc
    mechanism. See
    :class:`~sage.rpc.core.common.RemoteProcedureCaller` for a
    real-world example.

    INPUT:

    - ``instance`` -- the instance of the class that uses the
      ``@remote_callable`` decorator.

    OUTPUT:

    This iterater yields pairs consisting of the RPC name (string) and
    the RPC implementation (function or bound method).
    """
    seen = set()
    for name in dir(instance):
        holder = getattr(instance, name)
        if not isinstance(holder, remote_callable):
            continue
        if holder.name in seen:
            raise RuntimeError('RPC method name "{0}" already in use'.format(holder.name))
        else:
            seen.add(holder.name)
        if sys.version_info.major < 3:
            method = types.MethodType(holder.method, instance, instance.__class__)
        else:
            method = types.MethodType(holder.method, instance)
        setattr(instance, name, method)
        yield holder.name, method


class DecoratorTest(object):

    def __init__(self):
        """
        Test class for the :class:`remote_callable` decorator.

        The constructor shows how to iterate over the decorated
        methods. In real use, they would be registered in the RPC
        mechanism to be callable at the remote end. For purposes of
        this test, they are only printed.

        EXAMPLES::

            sage: from sage.rpc.core.decorator import DecoratorTest
            sage: test = DecoratorTest()
            Register foo to call _impl_foo
        """
        for name, method in remote_callable_iter(self):
            print('Register {0} to call {1}'.format(name, method.__name__))

    @remote_callable('foo')
    def _impl_foo(self):
        """
        The implementation of the RPC call.

        As a convention, the implementations of rpc methods start with
        ``'_impl_'``, which makes it clear that they are not supposed
        to be called locally. You are encouraged to follow this
        convention in your code as well.

        EXAMPLES::

            sage: from sage.rpc.core.decorator import DecoratorTest
            sage: test = DecoratorTest()
            Register foo to call _impl_foo
            sage: test._impl_foo()
            'Implementation'
        """
        return 'Implementation'

    
