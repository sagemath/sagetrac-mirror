r"""
The Cython bits of the citation management.

Cython functions to enable citation management, and to register citations.

AUTHOR:

- Martin Westerholt-Raum (martin@raum-brothers.eu)
"""

from sage.misc.citation_items.citation_item_dispatcher import get_citation_item

cdef bint sage_citation_enabled = False
cdef list sage_citation_records = list()

cpdef citation_add_record(record):
    r"""
    Add a citation record to the list of records citations will be
    stored in.

    INPUT:

    - ``record`` -- A citation record.

    EXAMPLES::

        sage: from sage.misc.cite import *
        sage: record = CitationRecord()
        sage: citation_add_record(record)
        sage: citation_enable()
        sage: cite("gmp")
        sage: record
        [GMP]
        sage: citation_disable(record)
    """
    global sage_citation_records
    sage_citation_records.append(record)

cpdef citation_remove_record(record):
    r"""
    Remove a citation record from the list of records citations will
    be stored in.

    INPUT:

    - ``record`` -- A citation record.

    EXAMPLES::

        sage: from sage.misc.cite import *
        sage: record = CitationRecord()
        sage: citation_enable(record)
        sage: citation_remove_record(record)
        sage: cite("gmp")
        sage: record
        []
        sage: citation_disable()
    """
    global sage_citation_records
    sage_citation_records.remove(record)

cpdef citation_enable(record = None):
    r"""
    Enable tracking of citations, adding ``record`` to the internal
    list of citation records.

    INPUT:

    - ``record`` -- A citation record or ``None``.

    EXAMPLES::

        sage: from sage.misc.cite import *
        sage: record = CitationRecord()
        sage: citation_enable(record)
        sage: cite("gmp")
        sage: citation_disable()
        sage: record
        [GMP]
    """
    global sage_citation_enabled

    if record is not None:
        citation_add_record(record)

    sage_citation_enabled = True

cpdef citation_disable(record = None):
    r"""
    Disable tracking of citations, removing ``record`` from the
    internal list of citation records.

    INPUT:

    - ``record`` -- A citation record or ``None``.

    EXAMPLES::

        sage: from sage.misc.cite import *
        sage: record = CitationRecord()
        sage: citation_enable(record)
        sage: cite("gmp")
        sage: citation_disable()
        sage: cite("maxima")
        sage: record
        [GMP]
    """
    global sage_citation_enabled
    sage_citation_enabled = False

    if record is not None:
        citation_remove_record(record)

cpdef cite(char* item_str):
    r"""
    Cite an item, using it's string.

    INPUT:

    - ``item_str`` -- A string.

    EXAMPLES::

        sage: from sage.misc.cite import *
        sage: record = CitationRecord()
        sage: citation_enable(record)
        sage: cite("gmp")
        sage: citation_disable(record)
        sage: record
        [GMP]
    """
    ## This is crucial to make use of branch prediction.  It results
    ## in vitually zero costs for calling this from cython.
    if not sage_citation_enabled:
        return
    
    cdef item = get_citation_item(item_str)
    for record in sage_citation_records:
        if item not in record:
            record.append(item)

