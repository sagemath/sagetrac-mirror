from sage.structure.element cimport Element, RingElement
from sage.rings.morphism cimport RingHomomorphism
from sage.rings.morphism cimport RingHomomorphism_im_gens


cdef class BoundedSeriesBaseringInjection(RingHomomorphism):
    cdef RingElement _an_element
    cdef object _new_constant_series_

cdef class BoundedSeriesRestriction(RingHomomorphism):
    cdef _series_class
    cdef _diff_log_radius

cdef class BoundedSeriesHomomorphism_im_gens(RingHomomorphism_im_gens):
    cdef _morphism
    cdef _image
    cdef _e
    cdef _function
    cdef _zeroes
    cdef _gain_precision
    cdef list _vertices
