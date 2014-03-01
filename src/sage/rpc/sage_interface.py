"""
Interface to Sage

This example shows how to write an interface to a remote Sage
session using the RPC framework.

EXAMPLES::

    sage: s = sage_remote.test_blocking()
    sage: s.eval('1+1')
    '2'
"""

import os

from sage.structure.sage_object import SageObject
from sage.rpc.sage_client import SageClient


class SageInterfaceBase(SageObject):

    def __init__(self, sage_client):
        """
        Base class for interface to Sage

        INPUT:

        - ``sage_client`` --
          :class:`sage.rpc.sage_client.SageClient`. RPC client for the
          sage compute service.

        EXAMPLES::

            sage: sage_remote.test_blocking()
            Blocking Sage interface
        """
        assert isinstance(sage_client, SageClient)
        self._client = sage_client
        self._client.interface = self
        self._current_label = 0
        self._replies = dict()

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::
        
            sage: s = sage_remote.test_blocking()
            sage: s._repr_()
            'Blocking Sage interface'
        """
        return 'Blocking Sage interface'

    def _new_label(self):
        """
        Return a new label

        OUTPUT:

        A unique JSON-serializable label.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: label = s._new_label()
            sage: isinstance(label, int)
            True
        """
        self._current_label += 1
        return self._current_label

    def _push(self, label, reply):
        """
        Temporarily save a reply
        
        INPUT:

        - ``label`` -- anything that identifies the result.

        - ``reply`` -- anything.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s._push('my_label', '123456')
            sage: s._pop('my_label')
            '123456'
        """
        self._replies[label] = reply

    def _pop(self, label):
        """
        Return a previously-pushed result.

        INPUT:

        - ``label`` -- anything that identifies a
          previouly-:meth:`_push` ed result..

        OUTPUT:

        The previously-pushed result.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s._push('my_label', '123456')
            sage: s._pop('my_label')
            '123456'
        """
        return self._replies.pop(label)

    def _assert_no_pending_replies(self):
        """
        Verify that there are no pending results pushed.

        This method triggers an assertion unless all previously-pushed
        results have been popped.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s._assert_no_pending_replies()

            sage: s._push('my_label', '123456')
            sage: s._assert_no_pending_replies()
            Traceback (most recent call last):
            ...
            AssertionError

            sage: s._pop('my_label')
            '123456'
            sage: s._assert_no_pending_replies()            
        """
        assert len(self._replies) == 0


class BlockingSageInterface(SageInterfaceBase):
    """
    Blocking Interface to Sage
    
    Blocking means that calls wait for the result of an evaluation
    request and return the relpy immediately. This is mostly useful
    for testing purposes. In real examples, you would be using the
    asynchronous interface to do something else until the result has
    been computed by the remote Sage session.
    """

    def _impl_sage_eval_stdin(self, label):
        assert(False)

    def _impl_sage_eval_stdout(self, stdout, label):
        self._stdout += str(stdout)

    def _impl_sage_eval_stderr(self, stderr, label):
        self._stdout += str(stderr)

    def _impl_sage_eval_result(self, cpu_time, wall_time, label):
        if self._stderr:
            result = self._stdout + self._stderr.rstrip()
        else:
            result = self._stdout.rstrip()
        self._push(label, result)

    def _impl_sage_eval_crash(self, label):
        self._push(label, 'remote session disconnected unexpectedly')

    def _impl_code_completion_finished(self, base, completions, label):
        self._push(label, [base, completions])

    def _impl_pong(self, start_time, label):
        import time
        elapsed_sec = time.time() - start_time
        self._push(label, elapsed_sec)
    
    def ping(self):
        """
        Ping the remote Sage session

        OUTPUT:

        Floating point number. The elapsed time in seconds.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s.ping()    # random output
            0.040821075439453125
        """
        self._assert_no_pending_replies()
        label = self._new_label()
        import time
        self._client.rpc.util.ping(time.time(), label)
        self._client.wait('util.pong')
        result = self._pop(label)
        self._assert_no_pending_replies()
        return result
        
    def eval(self, sage_command):
        """
        Evaluate a Sage command remotely

        INPUT:
        
        - ``sage_command`` -- string.

        OUTPUT:

        The output of the remote Sage computation as a string.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s.eval('pi.n(200)')
            '3.1415926535897932384626433832795028841971693993751058209749'
        """
        self._assert_no_pending_replies()
        self._stdout = self._stderr = ''
        label = self._new_label()
        self._client.sage_eval(sage_command, label)
        while True:
            self._client.wait()
            try:
                result = self._pop(label)
                break
            except KeyError:
                pass
        self._assert_no_pending_replies()
        return result

    def complete(self, code_string, cursor_position):
        """
        Return code completion suggestions

        INPUT:

        - ``code_string`` -- string.

        - ``cursor_position`` -- integer.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()
            sage: s.complete('pi.n(10', 2)
            [u'%pinfo', u'%pinfo2', u'pi', u'pickle_function', u'piecewise']
            sage: s.complete('pi.n(10', 3)
            [u'pi.N', u'pi.Order', u'pi.abs', u'pi.add', u'pi.add_to_both_sides', ...]
            sage: s.complete('pi.n(10', 4)
            [u'pi.n', u'pi.negation', u'pi.nintegral', u'pi.nintegrate', ...]
        """
        self._assert_no_pending_replies()
        label = self._new_label()
        self._client.code_complete(code_string, int(cursor_position), label)
        self._client.wait('code_completion.finished')
        base, completions = self._pop(label)
        self._assert_no_pending_replies()
        return completions
