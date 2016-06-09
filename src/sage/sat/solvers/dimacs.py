"""
SAT-Solvers via DIMACS Files

Sage supports calling SAT solvers using the popular DIMACS format. This module implements
infrastructure to make it easy to add new such interfaces and some example interfaces.

Currently, interfaces to **RSat** [RS]_ and **Glucose** [GL]_ are included by default.

.. note::

    Our SAT solver interfaces are 1-based, i.e., literals start at 1. This is consistent with the
    popular DIMACS format for SAT solving but not with Pythion's 0-based convention. However, this
    also allows to construct clauses using simple integers.

AUTHORS:

- Martin Albrecht (2012): first version

Classes and Methods
-------------------
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

import os, sys, subprocess, shlex

from sage.sat.solvers.satsolver import SatSolver
from sage.misc.all import tmp_filename, get_verbose
from time import sleep

class DIMACS(SatSolver):
    """
    Generic DIMACS Solver.

    .. note::

        Usually, users won't have to use this class directly but some
        class which inherits from this class.

    .. automethod:: __init__
    .. automethod:: __call__
    """

    command = ""

    def __init__(self, command=None, filename=None, verbosity=0, **kwds):
        """
        Construct a new generic DIMACS solver.

        INPUT:

        - ``command`` - a named format string with the command to
          run. The string must contain {input} and may contain
          {output} if the solvers writes the solution to an output
          file. For example "sat-solver {input}" is a valid
          command. If ``None`` then the class variable ``command`` is
          used. (default: ``None``)

        - ``filename`` - a filename to write clauses to in DIMACS
          format, must be writable. If ``None`` a temporary filename
          is chosen automatically. (default: ``None``)

        - ``verbosity`` - a verbosity level, where zero means silent
          and anything else means verbose output. (default: ``0``)

        - ``**kwds`` - accepted for compatibility with other solves,
          ignored.

        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS()
            DIMACS Solver: ''
        """
        if filename is None:
            filename = tmp_filename()

        self._headname = filename
        self._verbosity = verbosity

        if command is not None:
            self._command = command
        else:
            self._command = self.__class__.command

        self._tail  = open(tmp_filename(),'w')
        self._var = 0
        self._lit = 0

    def __repr__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS(command="iliketurtles {input}")
            DIMACS Solver: 'iliketurtles {input}'
        """
        return "DIMACS Solver: '%s'"%(self._command)

    def __del__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: d = DIMACS(command="iliketurtles {input}")
            sage: del d
        """
        if not self._tail.closed:
            self._tail.close()
        if os.path.exists(self._tail.name):
            os.unlink(self._tail.name)

    def var(self, decision=None):
        """
        Return a *new* variable.

        INPUT:

        - ``decision`` - accepted for compatibility with other solvers, ignored.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
        """
        self._var+= 1
        return self._var

    def nvars(self):
        """
        Return the number of variables.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
            sage: solver.var(decision=True)
            2
            sage: solver.nvars()
            2
        """
        return self._var

    def add_clause(self, lits):
        """
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` - a tuple of integers != 0

        .. note::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
            sage: solver.var(decision=True)
            2
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver
            DIMACS Solver: ''
        """
        l = []
        for lit in lits:
            lit = int(lit)
            while abs(lit) > self.nvars():
                self.var()
            l.append(str(lit))
        l.append("0\n")
        self._tail.write(" ".join(l) )
        self._lit += 1

    def write(self, filename=None):
        """
        Write DIMACS file.

        INPUT:

        - ``filename`` - if ``None`` default filename specified at initialization is used for
          writing to (default: ``None``)

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: _ = solver.write()
            sage: for line in open(fn).readlines():
            ....:     print(line)
            p cnf 3 1
            1 -2 3 0

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, -2 , 3) )
            sage: _ = solver.write(fn)
            sage: for line in open(fn).readlines():
            ....:      print(line)
            p cnf 3 1
            1 -2 3 0
        """
        headname = self._headname if filename is None else filename
        head = open(headname, "w")
        head.truncate(0)
        head.write("p cnf %d %d\n"%(self._var,self._lit))
        head.close()

        tail = self._tail
        tail.close()

        head = open(headname,"a")
        tail = open(self._tail.name,"r")
        head.write(tail.read())
        tail.close()
        head.close()

        self._tail = open(self._tail.name,"a")
        return headname

    def clauses(self, filename=None):
        """
        Return original clauses.

        INPUT:

        - ``filename`` - if not ``None`` clauses are written to ``filename`` in
          DIMACS format (default: ``None``)

        OUTPUT:

            If ``filename`` is ``None`` then a list of ``lits, is_xor, rhs``
            tuples is returned, where ``lits`` is a tuple of literals,
            ``is_xor`` is always ``False`` and ``rhs`` is always ``None``.

            If ``filename`` points to a writable file, then the list of original
            clauses is written to that file in DIMACS format.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, 2, 3) )
            sage: solver.clauses()
            [((1, 2, 3), False, None)]

            sage: solver.add_clause( (1, 2, -3) )
            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 3 2
            1 2 3 0
            1 2 -3 0
            <BLANKLINE>
        """
        if filename is not None:
            self.write(filename)
        else:
            tail = self._tail
            tail.close()
            tail = open(self._tail.name,"r")

            clauses = []
            for line in tail.readlines():
                if line.startswith("p") or line.startswith("c"):
                    continue
                clause = []
                for lit in line.split(" "):
                    lit = int(lit)
                    if lit == 0:
                        break
                    clause.append(lit)
                clauses.append( ( tuple(clause), False, None ) )
            tail.close()
            self._tail = open(self._tail.name, "a")
            return clauses

    @staticmethod
    def render_dimacs(clauses, filename, nlits):
        """
        Produce DIMACS file ``filename`` from ``clauses``.

        INPUT:

        - ``clauses`` - a list of clauses, either in simple format as a list of
          literals or in extended format for CryptoMiniSat: a tuple of literals,
          ``is_xor`` and ``rhs``.

        - ``filename`` - the file to write to

        - ``nlits -- the number of literals appearing in ``clauses``

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, 2, -3) )
            sage: DIMACS.render_dimacs(solver.clauses(), fn, solver.nvars())
            sage: print(open(fn).read())
            p cnf 3 1
            1 2 -3 0
            <BLANKLINE>

        This is equivalent to::

            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 3 1
            1 2 -3 0
            <BLANKLINE>

        This function also accepts a "simple" format::

            sage: DIMACS.render_dimacs([ (1,2), (1,2,-3) ], fn, 3)
            sage: print(open(fn).read())
            p cnf 3 2
            1 2 0
            1 2 -3 0
            <BLANKLINE>
        """
        fh = open(filename, "w")
        fh.write("p cnf %d %d\n"%(nlits,len(clauses)))
        for clause in clauses:
            if len(clause) == 3 and clause[1] in (True, False) and clause[2] in (True,False,None):
                lits, is_xor, rhs = clause
            else:
                lits, is_xor, rhs = clause, False, None

            if is_xor:
                closing = lits[-1] if rhs else -lits[-1]
                fh.write("x" + " ".join(map(str, lits[:-1])) + " %d 0\n"%closing)
            else:
                fh.write(" ".join(map(str, lits)) + " 0\n")
        fh.close()

    def __call__(self, assumptions=None):
        """
        Run 'command' and collect output.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        TESTS:

        This class is not meant to be called directly::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver()
            Traceback (most recent call last):
            ...
            ValueError: No SAT solver command selected.
        """
        if assumptions is not None:
            raise NotImplementedError("Assumptions are not supported for DIMACS based solvers.")

        self.write()

        output_filename = None
        self._output = []

        command = self._command.strip()

        if not command:
            raise ValueError("No SAT solver command selected.")


        if "{output}" in command:
            output_filename = tmp_filename()
        command = command.format(input=self._headname, output=output_filename)

        args = shlex.split(command)

        try:
            process = subprocess.Popen(args, stdout=subprocess.PIPE)
        except OSError:
            raise OSError("Could run '%s', perhaps you need to add your SAT solver to $PATH?"%(" ".join(args)))

        try:
            while process.poll() is None:
                for line in iter(process.stdout.readline,''):
                    if get_verbose() or self._verbosity:
                        print(line)
                        sys.stdout.flush()
                    self._output.append(line)
                sleep(0.1)
            if output_filename:
                self._output.extend(open(output_filename).readlines())
        except BaseException:
            process.kill()
            raise

class RSat(DIMACS):
    """
    An instance of the RSat solver.

    For information on RSat see: http://reasoning.cs.ucla.edu/rsat/
    """

    command = "rsat {input} -v -s"

    def __call__(self, assumptions=None):
        """
        Solve this instance.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLE::

           sage: from sage.sat.boolean_polynomials import solve as solve_sat
           sage: F,s = mq.SR(1,1,1,4,gf2=True,polybori=True).polynomial_system()
           sage: solve_sat(F, solver=sage.sat.solvers.RSat)  # optional - RSat
        """
        DIMACS.__call__(self)

        s = [None] + [False for _ in range(self.nvars())]
        for line in self._output:
            if line.startswith("c"):
                continue
            if line.startswith("s"):
                if "UNSAT" in line:
                    return False
            if line.startswith("v"):
                lits = map(int, line[2:-2].strip().split(" "))
                for e in lits:
                    s[abs(e)] = e>0
        return tuple(s)

class Glucose(DIMACS):
    """
    An instance of the Glucose solver.

    For information on Glucose see: http://www.lri.fr/~simon/?page=glucose
    """

    command = "glucose_static -verb=2 {input} {output}"

    def __call__(self, **kwds):
        """
        Solve this instance.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLE::

           sage: from sage.sat.boolean_polynomials import solve as solve_sat
           sage: F,s = mq.SR(1,1,1,4,gf2=True,polybori=True).polynomial_system()
           sage: solve_sat(F, solver=sage.sat.solvers.Glucose)  # optional - Glucose
        """
        DIMACS.__call__(self)

        for line in self._output:
            if line.startswith("c"):
                continue
            if line.startswith("s"):
                if "UNSAT" in line:
                    return False
            try:
                s = map(int, line[:-2].strip().split(" "))
                s = (None,) + tuple(e>0 for e in s)
                return s
            except ValueError:
                pass
        return False
