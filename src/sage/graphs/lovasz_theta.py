# -*- coding: utf-8 -*-
r"""
Lovász theta-function of graphs

AUTHORS:

- Dima Pasechnik (2015-06-30): Initial version

REFERENCE:

.. [Lovasz1979] László Lovász,
  "On the Shannon capacity of a graph",
  IEEE Trans. Inf. Th. 25(1979), 1-7.

Functions
---------
"""

def lovasz_theta(graph):
    r"""
    Return the value of Lovász theta-function of graph

    For a graph `G` this function is denoted by `\theta(G)`, and it can be
    computed in polynomial time. Mathematically, its most important property is the following:

    .. MATH::

        \alpha(G)\leq\theta(G)\leq\chi(\overline{G})

    with `\alpha(G)` and `\chi(\overline{G})` being, respectively, the maximum
    size of an :meth:`independent set <sage.graphs.graph.Graph.independent_set>`
    set of `G` and the :meth:`chromatic number
    <sage.graphs.graph.Graph.chromatic_number>` of the :meth:`complement
    <sage.graphs.generic_graph.GenericGraph.complement>` `\overline{G}` of `G`.

    For more information, see the :wikipedia:`Lovász_number`.

    .. NOTE::

        - Implemented for undirected graphs only. Use to_undirected to convert a
          digraph to an undirected graph.

        - This function requires the optional package ``csdp``, which you can
          install with with ``sage -i csdp``.

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: C.lovasz_theta()                             # optional csdp
          4.0
          sage: graphs.CycleGraph(5).lovasz_theta()          # optional csdp
          2.236068

    TEST::

        sage: g = Graph()
        sage: g.lovasz_theta() # indirect doctest
        0
    """
    n = graph.order()
    if n == 0:
        return 0

    from networkx import write_edgelist
    from sage.misc.temporary_file import tmp_filename
    import os, subprocess

    CSDP().require()

    g = graph.relabel(inplace=False, perm=range(1,n+1)).networkx_graph()
    tf_name = tmp_filename()
    tf = open(tf_name, 'wb')
    tf.write(str(n)+'\n'+str(g.number_of_edges())+'\n')
    write_edgelist(g, tf, data=False)
    tf.close()
    lines = subprocess.check_output(['theta', tf_name])
    return float(lines.split()[-1])

from sage.misc.feature import Executable
class CSDP(Executable):
    r"""
    A class:`sage.misc.feature.Feature` which checks for the ``theta`` binary
    of CSDP.

    EXAMPLES::

        sage: from sage.graphs.lovasz_theta import CSDP
        sage: CSDP().is_present() # optional: csdp
        True

    """
    def __init__(self):
        Executable.__init__(self, name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp")

    def is_functional(self):
        r"""
        Check whether ``theta`` works on a trivial example.

        EXAMPLES::

            sage: from sage.graphs.lovasz_theta import CSDP
            sage: CSDP().is_functional() # optional: csdp
            True

        """
        from sage.misc.temporary_file import tmp_filename
        from sage.misc.misc import verbose
        import os, subprocess
        tf_name = tmp_filename()
        with open(tf_name, 'wb') as tf:
            tf.write("2\n1\n1 1")
        try:
            lines = subprocess.check_output(['theta', tf_name])
        except subprocess.CalledProcessError as e:
            verbose("Call to theta failed with exit code %s and output:\n%s"%(e.returncode, lines), level=Feature.VERBOSE_LEVEL)
            return False
        result = lines.strip().split('\n')[-1]
        import re
        match = re.match("^The Lovasz Theta Number is (.*)$", result)
        if match is None:
            verbose("Last line of theta output did not match the expected format: `%s`"%(result,), level=Feature.VERBOSE_LEVEL)
            return False
        return True
