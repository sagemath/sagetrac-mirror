r"""
Interface to Msieve
\url{http://sourceforge.net/projects/msieve/}
"""
import os

from sage.rings.integer import Integer
from sage.misc.all import tmp_dir
from sage.env import SAGE_LOCAL

binary = os.path.join(SAGE_LOCAL,"bin","msieve")

_tmp_dir_existing = False
_tmp_dir_name = ""

def tmpdir():
    r"""
    Create temporary directory for msieve.
    """
    global _tmp_dir_existing
    global _tmp_dir_name
    if _tmp_dir_existing:
        return
    _tmp_dir_name = tmp_dir('msieve')
    _tmp_dir_existing = True

def msieve(number, block=True, time=False, verbose=False, clear=True, **kwds):
    r"""
    Return factorization of ``number``.

    INPUT:

        ``number`` -- an integer
            numbers starting with '0' are treated as octal,
            numbers starting with '0x' are treated as hexadecimal
        ``block`` -- (default: True) if True, you must wait until the
            sieve computation is complete before using Sage further.
            If False, Sage will run while the sieve computation
            runs in parallel.  If m is the returned object, use
            m.quit() to terminate a running factorization.
        ``time`` -- (default: False) if True, also used time is returned
        ``verbose`` -- (default: False) if True, print out verbose
            logging information about what happened during
            the Sieve run (for non-blocking Sieve, verbose information
            is always available via the log() method.)
        ``clear`` -- (default: True) if False, msieve uses previous results
            for following function calls

        s -- name   save intermediate results to <name>
                    instead of the default msieve.dat
        l -- name   append log information to <name>
                    instead of the default msieve.log
        i -- name   read one or more integers to factor from
                    <name> (default worktodo.ini) instead of
                    from the command line
        m --        manual mode: enter numbers via standard input
        q --        quiet: do not generate any log information,
                    only print any factors found
        d -- min    deadline: if still sieving after <min>
                    minutes, shut down gracefully (default off)
        r -- num    stop after finding <num> relations
        v --        verbose: write log information to screen
                    as well as to logfile
        t -- num    use at most <num> threads

        elliptic curve options:
        e --            perform 'deep' ECM, seek factors > 15 digits
        quadratic sieve options:
        c --            client: only perform sieving

        number field sieve options:
        n --            use the number field sieve (97+ digits only;
                        performs all NFS tasks in order)
        nf -- name      read from / write to NFS factor base file
                        <name> instead of the default msieve.fb
        np -- [X,Y]     perform only NFS polynomial selection; if
                        specified, cover the range from X to Y
                        percent of the total polynomial search
                        space, inclusive (1<=X<Y<=100)
        ns -- [X,Y]     perform only NFS sieving; if specified,
                        handle sieve lines X to Y inclusive
        nc --           perform only NFS combining (all phases)
        nc1 -- [X,Y]    perform only NFS filtering. Filtering will
                        only track rational ideals >= X and algebraic
                        ideals >= Y (determined automatically if 0 or
                        unspecified)
        nc2 --          perform only NFS linear algebra
        ncr --          perform only NFS linear algebra, restarting
                        from a previous checkpoint
        nc3 -- [X,Y]    perform only NFS square root (compute
                        dependency numbers X through Y, 1<=X,Y<=64)

    OUTPUT:

        ``list`` -- a list of the prime factors of n found
        ``str`` -- the time in cpu seconds that the computation took.
               (If time is False, this is always an empty string.)

    EXAMPLES::

        sage: msieve(3444)                      # optional -- msieve
        ([2, 2, 3, 7, 41], '')
        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1)) # optional -- msieve
        sage: msieve(n)                         # optional -- msieve
        ([10000000000000000051, 100000000000000000039], '')
        sage: msieve(n, time=True)              # random optional -- msieve

        sieving in progress (press Ctrl-C to pause)
        879 relations (348 full + 531 combined from 3127 partial), need 602
        sieving complete, commencing postprocessing
        ([10000000000000000051, 100000000000000000039], 'elapsed time 00:00:00')

        sage: msieve(3444, verbose=True)        # random optional -- msieve
        <BLANKLINE>
        <BLANKLINE>
        Msieve v. 1.49
        Tue Nov 1 18:09:48 2011
        random seeds: f8d07204 129e4eb9
        factoring 3444 (4 digits)
        p1 factor: 2
        p1 factor: 2
        p1 factor: 3
        p1 factor: 7
        p2 factor: 41
        elapsed time 00:00:00
        <BLANKLINE>
        ([2, 2, 3, 7, 41], '')
        sage: msieve(n, n='', nc3='[5,32]')     # optional -- msieve
        ([10000000000000000051, 100000000000000000039], '')

    AUTHORS:
        Jason Papadopoulos wrote msieve.

    Interface written by
    - Jeffrey Blakeslee (2009-02-18): initial version
    - Andre Apitzsch (2010-10-27): documentation and clear option
    """

    if not os.path.exists(binary):
        raise NotImplementedError("You might consider installing the optional package msieve with the following command: ``sage -i msieve-1.49.p0``")
        return
    if block:
        return msieve_block(number, time, verbose, clear, kwds)
    else:
        return msieve_nonblock(number, time, verbose, clear, kwds)

def msieve_block(n, time, verbose, clear, kwds):
    r"""
    Return factorization of n.
    """
    tmpdir()
    if clear and os.path.exists(_tmp_dir_name + '/msieve.dat'):
        os.remove(_tmp_dir_name + '/msieve.dat')
    options = ' '.join(['-%s %s'%(x,v) for x, v in kwds.iteritems()])
    if verbose or time:
        out = os.popen('msieve -m -v -s %s -l %s %s %s'%(_tmp_dir_name +
            '/msieve.dat', _tmp_dir_name + '/msieve.log', options, n)).read()
    else:
        out = os.popen('msieve -m -q -s %s -l %s %s %s'%(_tmp_dir_name +
            '/msieve.dat', _tmp_dir_name + '/msieve.log', options, n)).read()
    factorlist, t, verb = data_to_list(out, n, time=time, verbose=verbose)

    if verbose:
        print verb
    if time:
        return factorlist, t
    return factorlist, ''

def data_to_list(out, n, time, verbose):
    r"""
    Select factors from msieve output.
    """
    i = out.find('failed to reallocate')
    if i != -1:
        raise MemoryError, out
    out = out.replace('error:','error ')
    if verbose or time:
        l = out.split('factor:')
    else:
        l = out.split(':')
    result =[]
    Z = Integer
    for i in range(1,len(l)):
        result.append(Z(l[i].split("\n")[0]))
    t = l[len(l)-1].split("\n")[1].replace("\r","")
    return result, t, out

import pexpect
import cleaner
class msieve_nonblock:
    r"""
    A non-blocking version of msieve.

    The msieve starts running when you create the object, but you can
    still use Sage in parallel:

    EXAMPLES::

        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: m = msieve(n, block=False, time=True) # optional -- msieve
        sage: m                                     # random optional -- msieve
        Proper factors so far: []
        sage: m                                     # random optional -- msieve
        Proper factors so far: [10000000000000000051, 100000000000000000039]
        sage: m                                     # random optional -- msieve
        ([10000000000000000051, 100000000000000000039], 'elapsed time 00:00:03')
        sage: m.list()                              # random optional -- msieve
        [10000000000000000051, 100000000000000000039]
        sage: m.time()                              # random optional -- msieve
        'elapsed time 00:00:03'

        sage: m = msieve(next_prime(10^20)*next_prime(10^21), block=False)  # optional -- msieve
        sage: m                                     # random optional -- msieve
        Proper factors so far: [100000000000000000039, 1000000000000000000117]
        sage: m                                     # random optional -- msieve
        [100000000000000000039, 1000000000000000000117]

        sage: k = 19; n = next_prime(10^k)*next_prime(10^(k+1))
        sage: m = msieve(n, block=False, time=True) # optional -- msieve requires time command
        sage: m                                     # random optional -- msieve
        Proper factors so far: []
        sage: m                                     # random optional -- msieve
        ([10000000000000000051, 100000000000000000039], '0.21')
        sage: m.list()                              # random optional -- msieve
        [10000000000000000051, 100000000000000000039]
        sage: m.time()                              # random optional -- msieve requires time command
        '0.21'

        sage: m = msieve(next_prime(10^20)*next_prime(10^21), block=False)  # optional -- msieve
        sage: m                                     # random optional -- msieve
        Proper factors so far: [100000000000000000039, 1000000000000000000117]
        sage: m                                     # random optional -- msieve
        [100000000000000000039, 1000000000000000000117]
    """
    def __init__(self, n, time, verbose, clear, kwds):
        self._n = n
        tmpdir()
        if clear and os.path.exists(_tmp_dir_name + '/msieve.dat'):
            os.remove(_tmp_dir_name + '/msieve.dat')
        options = ' '.join(['-%s %s'%(x,v) for x, v in kwds.iteritems()])
        if verbose or time:
            cmd = 'msieve -m -v -s %s -l %s %s %s'%(_tmp_dir_name +
                  '/msieve.dat', _tmp_dir_name + '/msieve.log' ,options, n)
        else:
            cmd = 'msieve -m -q -s %s -l %s %s %s'%(_tmp_dir_name +
                  '/msieve.dat', _tmp_dir_name + '/msieve.log' ,options, n)
        self._p = pexpect.spawn(cmd)
        cleaner.cleaner(self._p.pid, 'msieve')
        self._done = False
        self._out = ''
        self._time = ''
        self._do_time = time
        self._do_verbose = verbose

    def n(self):
        r"""
        Return the integer that is being factored.
        """
        return self._n

    def pid(self):
        r"""
        Return the PIN id of the msieve process (actually
        of the time process that spawns the sieve process).
        """
        return self._p.pid

    def done(self):
        r"""
        Return True if the sieve process has completed.
        """
        return self._done

    def __repr__(self):
        r"""
        Return a text representation of self.
        """
        if self._done:
            if hasattr(self, '_killed') and self._killed:
                return "Factorization was terminated early."
            factorlist, t, verb = data_to_list(self._get(), self._n, self._do_time,self._do_verbose)
            if self._do_verbose:
                print verb
            if self._do_time:
                return str( (factorlist,t) )
            return str(factorlist)
        else:
            return 'Proper factors so far: %s'%self.list()

    def cputime(self):
        r"""
        Return the time in seconds (as a string) that it took to
        factor n, or return '?' if the factorization has not
        completed or the time is unknown.
        """
        if not self._do_time:
            raise ValueError, "you have to start the seive with the option time=True in order to get timing information"
        try:
            return data_to_list(self._get(), self._n, self._do_time, self._do_verbose)[1]
        except IndexError:
            return '?'
    time = cputime

    def log(self):
        r"""
        Return all output of running the sieve so far.
        """
        return self._get()

    def __getitem__(self, i):
        r"""
        Return the i-th factor (in sorted order) found so far.
        """
        return self.list()[i]

    def __len__(self):
        r"""
        Return the number of factors found so far.  If m is the
        Sieve object, type len(m) to see the number of factors.
        """
        return len(self.list())

    def list(self):
        r"""
        Return the a list of the factors found so far, as Sage
        integers.
        """
        try:
            return data_to_list(self._get(), self._n, self._do_time, self._do_verbose)[0]
        except IndexError:
            return []

    def quit(self):
        r"""
        Terminate the Msieve process, in case you want
        to give up on computing this factorization.

        EXAMPLES::

            sage: n = next_prime(2^110)*next_prime(2^100)
            sage: ms = msieve(n, block=False)   # optional -- msieve
            sage: ms                            # optional -- msieve
            Proper factors so far: []
            sage: ms.quit()                     # optional -- msieve
            sage: ms                            # optional -- msieve
            Factorization was terminated early.
        """
        pid = self.pid()
        os.killpg(int(pid),9)
        #self._p.close()
        self._killed = True
        self._done = True

    def _get(self, timeout=0.1):
        r"""
        Used internally to get information about what has been
        computed so far.
        """
        if self._done:
            return self._out
        e = self._p
        try:
            e.expect('xxx', timeout=timeout)
        except pexpect.TIMEOUT, msg:
            pass
        except pexpect.EOF, msg:
            pass
            self._done = True
            self._p.close()
        self._out += e.before
        return self._out

