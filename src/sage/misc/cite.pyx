r"""
The Cython bits of the citation management.

Cython functions to enable citation management, and to register citations.

AUTHOR:

- Martin Westerholt-Raum (martin@raum-brothers.eu)
"""

from sage.misc.citation_items.citation_item_dispatcher import get_citation_item

cdef bint sage_citation_enabled = False
cdef list sage_citation_records = list()

cpdef citation_enable(record = None):
    if record is not None:
        sage_citation_records.append(record)
    sage_citation_enabled = True

cpdef citation_disable(record = None):
    sage_citation_enabled = False
    if record is not None:
        sage_citation_records.remove(record)

cpdef cite(char* item_str):
    if not sage_citation_enabled:
        return
    
    cdef item = get_citation_item(item_str)
    for record in sage_citation_records:
        record.append(item)
