r"""
ReST index of functions

This module contains a function that generates a ReST index table of functions
for use in doc-strings.

{INDEX_OF_FUNCTIONS}

"""
def gen_rest_table_index(list_of_entries, sort=True, only_local=True):
    r"""
    Return a ReST table describing a list of functions.

    The list of functions can either be given explicitly, or implicitly as the
    functions/methods of a module or class.

    In the case of a class, only static non-inherited methods are listed.

    INPUT:

    - ``list_of_entries`` -- a list of functions, a module or a class. If given
      a list of functions, the generated table will consist of these. If given a
      module, all functions implemented in that module will be listed, except
      deprecated ones or those starting with an underscore. If a class is given,
      all non-inherited static methods of that class will be listed, except
      deprecated ones or those starting with an underscore.

    - ``sort`` (boolean; ``True``) -- whether to sort the list of methods
      lexicographically.

    - ``only_local`` (boolean; ``True``) -- if ``list_of_entries`` is a module,
       ``only_local = True`` means that imported functions will be filtered out.
       This can be useful to disable for making indexes of e.g. catalog modules
       such as ``sage.coding.codes_catalog``.

    EXAMPLE::

        sage: from sage.misc.rest_index_of_methods import gen_rest_table_index
        sage: print gen_rest_table_index([graphs.PetersenGraph])
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: |
        <BLANKLINE>
           :func:`~sage.graphs.generators.smallgraphs.PetersenGraph` | The Petersen Graph is a named graph that consists of 10 vertices...

    The table of a module::

        sage: print gen_rest_table_index(sage.misc.rest_index_of_methods)
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: |
        <BLANKLINE>
           :func:`~sage.misc.rest_index_of_methods.gen_rest_table_index` | Return a ReST table describing a list of functions...

    The table of a class::

        sage: print gen_rest_table_index(Graph)
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: |
        ...
           :meth:`~sage.graphs.graph.Graph.sparse6_string` | Returns the sparse6 representation of the graph as an ASCII string.
        ...

    TESTS:

    The inherited methods do not show up::

        sage: gen_rest_table_index(sage.combinat.posets.lattices.FiniteLatticePoset).count('\n') < 50
        True
        sage: from sage.graphs.generic_graph import GenericGraph
        sage: A = gen_rest_table_index(Graph).count('\n')
        sage: B = gen_rest_table_index(GenericGraph).count('\n')
        sage: A < B
        True

    When only_local is False, don't include gen_rest_table_index itself::

        sage: print gen_rest_table_index(sage.misc.rest_index_of_methods, only_local=True)
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: |
        <BLANKLINE>
           :func:`~sage.misc.rest_index_of_methods.gen_rest_table_index` | Return a ReST table describing a list of functions...
        sage: print gen_rest_table_index(sage.misc.rest_index_of_methods, only_local=False)
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: |
        <BLANKLINE>
        <BLANKLINE>
        
    """
    import inspect

        
        
    # If input is a class/module, we list all its non-private and methods/functions
    if (inspect.isclass(list_of_entries) or
        inspect.ismodule(list_of_entries)):
        root = list_of_entries
        def local_filter(f,name):
            module = inspect.getmodule(f)
            if only_local:
                return inspect.getmodule(root) == module
            else:
                return inspect.isclass(list_of_entries) or module.__name__ != "sage.misc.rest_index_of_methods"
        list_of_entries = [getattr(root,name) for name,f in root.__dict__.items() if
                           (not name.startswith('_')     and # private functions
                            not hasattr(f,'trac_number') and # deprecated functions
                            not inspect.isclass(f)       and # classes
                            local_filter(f,name)             # possibly filter imported functions
                            )]

    assert isinstance(list_of_entries,list)

    s = (".. csv-table::\n"
         "   :class: contentstable\n"
         "   :widths: 30, 70\n"
         "   :delim: |\n\n")

    if sort:
        list_of_entries.sort(key=lambda x:getattr(x,'__name__',''))

    for e in list_of_entries:

        if inspect.ismethod(e):
            link = ":meth:`~"+str(e.im_class.__module__)+"."+str(e.im_class.__name__)+"."+e.__name__+"`"
        elif inspect.isfunction(e):
            link = ":func:`~"+str(e.__module__)+"."+str(e.__name__)+"`"
        else:
            continue

        # Descriptions of the method/function
        if e.__doc__:
            desc = e.__doc__.splitlines()
            desc = desc[0] if desc[0] else desc[1]
        else:
            desc = "NO DOCSTRING"

        s += "   {} | {}\n".format(link,desc.lstrip())

    return s+'\n'

__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index([gen_rest_table_index]))
