r"""
Interface to Ruby

This is an expect interface to Ruby.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import cPickle
import os
import time
from expect import Expect, ExpectElement, FunctionElement
from sage.structure.sage_object import dumps, loads, load


class Ruby(Expect):
    r"""
    Expect interface to the Ruby interpreter.

    INPUT:

    - server -- (optional); if specified runs SAGE on a remote
                  machine with address.  You must have ssh keys setup
                  so you can login to the remote machine by typing
                  "ssh remote_machine" and no password, call
                  _install_hints_ssh() for hints on how to do that.

    EXAMPLES:

    We create an interface to Ruby.

        sage: from sage.interfaces.ruby import Ruby
        sage: ru = Ruby()

    Create the element 2 in our new copy of \sage, and cubeit. ::

        sage: a = ru(2)
        sage: a^3
        8

    Create a vector space of dimension $4$, and compute its generators::

        sage: V = ru('QQ^4')
        sage: V.gens()
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))

    Note that V is a not a vector space, it's a wrapper around an object
    (which happens to be a vector space), in another running instance
    of \sage. ::

        sage: type(V)
        <class 'sage.interfaces.ruby.RubyElement'>
        sage: V.parent()
        Ruby
        sage: g = V.0;  g
        (1, 0, 0, 0)
        sage: g.parent()
        Ruby

    We can still get the actual parent by using the name attribute of g,
    which is the variable name of the object in the child process. ::

        sage: ru('%s.parent()'%g.name())
        Vector space of dimension 4 over Rational Field

    Note that the memory space is completely different. ::

        sage: x = 10
        sage: ru('x = 5')
        5
        sage: x
        10
        sage: ru('x')
        5
    """
    def __init__(self, logfile=None,
                       python=False,
                       init_code=None,
                       server=None,
                       server_tmpdir=None,
                       remote_cleaner=True,
                       **kwds):
        command = "irb --noreadline"
        prompt = ">> "

        Expect.__init__(self,
                        name='ruby',
                        prompt=prompt,
                        command=command,
                        restart_on_ctrlc=False,
                        logfile=logfile,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        remote_cleaner=remote_cleaner,
                        **kwds
                        )

    def cputime(self, t=None):
        """
        Return cputime since this Ruby subprocess was started.

        EXAMPLES::

            sage: from sage.interfaces.ruby import Ruby
            sage: ruby = Ruby()
            sage: ruby.cputime()     # random output
            1.3530439999999999
            sage: ruby('factor(2^157-1)')
            852133201 * 60726444167 * 1654058017289 * 2134387368610417
            sage: ruby.cputime()     # random output
            1.6462939999999999
        """
        s = self.eval('cputime(%s)' % t)
        i = s.rfind('m')
        if i != -1:
            s = s[i + 1:-1]
        return float(s)

    def trait_names(self):
        return eval(self.eval('globals().keys()'))

    def quit(self, verbose=False):
        import signal
        if not self._expect is None:
            pid = self._expect.pid
            if verbose:
                if self.is_remote():
                    print "Exiting spawned %s process (local pid=%s, running on %s)" % (self, pid, self._server)
                else:
                    print "Exiting spawned %s process (pid=%s)." % (self, pid)
            try:
                for i in range(10):   
                    # multiple times, since clears out junk injected
                    # with ._get, etc.
                    self._expect.sendline(chr(3))  # send ctrl-c
                    self._expect.sendline('quit_sage(verbose=%s)' % verbose)
                    self._so_far(wait=0.2)

                os.killpg(pid, 9)
                os.kill(pid, 9)

            except (RuntimeError, OSError), msg:
                pass

            try:
                os.killpg(pid, 9)
                os.kill(pid, 9)
            except OSError:
                pass

            try:
                self._expect.close(signal.SIGQUIT)
            except Exception:
                pass
            self._expect = None

    def _remote_tmpfile(self):
        """
        """
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = eval(self.eval('import sage.interfaces.expect as e; e.tmp'))
            return self.__remote_tmpfile

    def __call__(self, x):
        """
        """
        if isinstance(x, RubyElement) and x.parent() is self:
            return x
        if isinstance(x, str):
            return RubyElement(self, x)

        if self.is_local():
            open(self._local_tmpfile(), 'w').write(cPickle.dumps(x, 2))
            return RubyElement(self, 'cPickle.load(open("%s"))' % self._local_tmpfile())
        else:
            open(self._local_tmpfile(),'w').write(dumps(x))   # my dumps is compressed by default
            self._send_tmpfile_to_server()
            return RubyElement(self, 'loads(open("%s").read())' % self._remote_tmpfile())

    def __reduce__(self):
        return reduce_load_Ruby, tuple([])

    def _quit_string(self):
        return 'from sage.misc.misc import delete_tmpfiles; delete_tmpfiles()'

    def eval(self, line, strip=True):
        """
        Send the code x to a second instance of the \sage interpreter and
        return the output as a string.

        This allows you to run two completely independent copies of \sage
        at the same time in a unified way.

        INPUT:

        - line -- input line of code
        - strip -- ignored
        """
        return Expect.eval(self, line)

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s' % (var, value)
        out = self.eval(cmd)
        if 'Traceback' in out:
            raise TypeError("Error executing code in SAGE\nCODE:\n\t%s\nSAGE ERROR:\n\t%s" % (cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.
        """
        s = self.eval('puts %s' % var)
        i = s.rfind('=>')
        if i != -1:
            return s[:i]
        return s

    #def clear(self, var):
    #    """
    #    Clear the variable named var.
    #    """
    #    self.eval('del %s'%var)

    def _contains(self, v1, v2):
        return self.eval('%s in %s' % (v1, v2))

    def _is_true_string(self, t):
        return t == "True"

    def console(self):
        ruby_console()

    def version(self):
        """
        """
        return eval(ruby_version())

    def _object_class(self):
        return RubyElement

    def new(self, x):
        return RubyElement(self, x)


class RubyElement(ExpectElement):
    def __getattr__(self, attrname):
        self._check_valid()
        return RubyFunction(self, attrname)

    def _sage_(self):
        """
        Return local copy of self.

        EXAMPLE::

            sage: from sage.interfaces.ruby import Ruby
            sage: ruby = Ruby()
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == ruby(F)._sage_()
            True
        """
        P = self.parent()
        if P.is_remote():
            P.eval('save(%s, "%s")' % (self.name(), P._remote_tmpfile()))
            P._get_tmpfile_from_server(self)
            return load(P._local_tmp_file())
        else:
            P.eval('save(%s, "%s")' % (self.name(), P._local_tmpfile()))
            return load(P._local_tmpfile())


class RubyFunction(FunctionElement):
    def __call__(self, *args, **kwds):
        P = self._obj.parent()
        args = [P(x) for x in args]
        args = ','.join([x.name() for x in args])
        kwds = ",".join(["%s=%s"%(k,P(v).name()) for k,v in kwds.iteritems()])
        if args != "" and kwds != "":
            callstr = '%s.%s(%s,%s)'%(self._obj._name, self._name, args, kwds)
        elif kwds != "":
            callstr = '%s.%s(%s)'%(self._obj._name, self._name, kwds)
        elif args != "":
            callstr = '%s.%s(%s)'%(self._obj._name, self._name, args)
        else:
            callstr = '%s.%s()'%(self._obj._name, self._name)
        z = RubyElement(P, callstr)

        return z

    def __repr__(self):
        return str(self.eval('%s.%s'%(self._obj._name, self._name)))


ruby = Ruby()


def reduce_load_ruby():
    return ruby


def reduce_load_element(s):
    return ruby('loads({})'.format(s))


def ruby_console():
    os.system('sage')


def ruby_version():
    return ruby('version()')
