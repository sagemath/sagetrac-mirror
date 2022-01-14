# -*- coding: utf-8 -*-
r"""
Cubic Hecke Database

This module contains the class :class:`CubicHeckeDataBase` which serves as an
interface to Ivan Marin's data files with respect to the cubic Hecke algebras.
Furthermore, it contains the class :class:`CubicHeckeFileCache` which enables
:class:`CubicHeckeAlgebra` to keep intermediate results of calculations in the
file system.

AUTHORS:

- Sebastian Oehms May 2020: initial version
"""


##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################


import os
from enum import Enum

from sage.structure.sage_object import SageObject
from sage.misc.persist import _base_dumps, save, load
from sage.misc.temporary_file import atomic_write
from sage.misc.verbose import verbose
from sage.env import SAGE_SHARE, SAGE_ROOT
from sage.matrix.constructor import matrix, Matrix  # uppercase version used in Marin's file `MatricesRegH4.maple`
from sage.rings.integer_ring import ZZ
from sage.algebras.hecke_algebras.base_rings_of_definition.cubic_hecke_base_ring import CubicHeckeRingOfDefinition





#------------------------------------------------------------------------------
# functions to convert matrices and ring elements to and from flat python
# dictionaries in order to save matrices avoiding compatibility problems with
# older or newer sage versions and to save disc space
#------------------------------------------------------------------------------
def simplify(elem):
    r"""
    Convert an element to a python dictionary recursively using its
    :meth:`_reconstruction_data`.

    INPUT:

    -- ``elem`` - element to be converted into python dictionary

    OUTPUT:

    A python dictionary from which ``elem`` can be reconstructed via
    element construction. The values of the dictionary may be
    dictionaries again.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import simplify
        sage: L.<c>=LaurentPolynomialRing(ZZ, 'c')
        sage: P.<a,b> = L['a,b']; P
        Multivariate Polynomial Ring in a, b
          over Univariate Laurent Polynomial Ring in c over Integer Ring
        sage: elem = 5*b-3*a*~c; elem
        (-3*c^-1)*a + 5*b
        sage: simplify(elem)
        {(0, 1): {0: 5}, (1, 0): {-1: -3}}
        sage: mat = matrix(P, [[2*a, -3], [c, 4*b*~c]]); mat
        [       2*a         -3]
        [         c (4*c^-1)*b]
        sage: simplify(mat)
        {(0, 0): {(1, 0): {0: 2}},
        (0, 1): {(0, 0): {0: -3}},
        (1, 0): {(0, 0): {1: 1}},
        (1, 1): {(0, 1): {-1: 4}}}
    """
    return elem._reconstruction_data()


class CubicHeckeDataSection(Enum):
    r"""
    Enum for the different sections of the database. The following choices are
    possible:

    - ``basis``  -- list of basis elements
    - ``reg_left_reprs``  -- data for the left regular representation
    - ``reg_right_reprs``  -- data for the right regular representation
    - ``irr_reprs`` -- data for the split irreducible representations

    Examples::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db.section
        <enum 'CubicHeckeDataSection'>
    """
    basis         = 'basis'
    regular_left  = 'regular_left'
    regular_right = 'regular_right'
    split_irred   = 'split_irred'


#-------------------------------------------------------------------------------
# Class to supply data for the basis and matrix representation for the cubic
# Hecke algebra
#-------------------------------------------------------------------------------
class CubicHeckeDataBase(SageObject):
    r"""
    Database interface needed for :class:`CubicHeckeAlgebra`.

    The original data are obtained from Ivan Marin's web-page (URL see the
    example below). In order to have these data installed during the build
    process as a sage-package they are converted as python files into a tarball.
    This tarball has been created using the method :meth:`create_spkg_tarball`.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db._feature
        Feature('database_cubic_hecke')
    """

    section = CubicHeckeDataSection

    def __init__(self):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db._data_library
            {}
        """
        from sage.features.databases import DatabaseCubicHecke
        self._feature   = DatabaseCubicHecke()
        self._data_library = {}
        self._demo = None

    def version(self):
        r"""
        Return the current version.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.version() > '2022.1.1'
            True
        """
        self._feature.require()
        from database_cubic_hecke import version
        return version()

    def demo_version(self):
        r"""
        Return whether the KnotInfo databases are installed completely or
        just the demo version is used.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.demo_version()       # optional - database_knotinfo
            False
        """
        if self._demo is None:
            if self._feature.is_present():
                self._demo = False
            else:
                self._demo = True
                self._data_library = demo_library
        return self._demo

    # --------------------------------------------------------------------------
    # read from an sobj-file obtained from Ivan Marin's database
    # --------------------------------------------------------------------------
    def read(self, section, translation_dict=None, nstrands=4):
        r"""
        Access various static data library.

        INPUT:

        ``section`` -- instance of enum :class:`CubicHeckeDataSection`
          to select the data to be read in

        OUTPUT:

        A dictionary containing the data corresponding to the section.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: basis = cha_db.read(cha_db.section.basis, nstrands=3)
            sage: len(basis)
            24
        """
        if not isinstance(section, CubicHeckeDataSection):
            raise TypeError('section must be an instance of enum %s' %(CubicHeckeDataBase.section))

        data_lib = self._data_library

        nstrands = int(nstrands)
        if (section, nstrands) in data_lib.keys():
            return data_lib[(section, nstrands)]

        verbose('loading data library %s for %s strands ...' %(section.value, nstrands))

        if self.demo_version():
            if nstrands >= 4:
                self._feature.require()
        else:
            from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import GenSign
            from database_cubic_hecke import read_basis, read_irr, read_reg
            if section == CubicHeckeDataSection.basis:
                data_lib[(section, nstrands)] = read_basis(num_strands=nstrands)
            elif section == CubicHeckeDataSection.split_irred:
                dim_list, repr_list, repr_list_inv = read_irr(translation_dict=translation_dict, num_strands=nstrands)
                data_lib[(section, nstrands)] = {GenSign.pos:repr_list, GenSign.neg:repr_list_inv}
            else:
                right = False
                if section == CubicHeckeDataSection.regular_right:
                    right = True
                dim_list, repr_list, repr_list_inv = read_reg(translation_dict=translation_dict, right=right, num_strands=nstrands)
                data_lib[(section, nstrands)] = {GenSign.pos:repr_list, GenSign.neg:repr_list_inv}

        verbose('... finished!')

        return data_lib[(section,nstrands)]


    # --------------------------------------------------------------------------
    # matrix_reprs_from_file_cache_
    # --------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, gen_ind, nstrands, ring_of_definition):
        r"""
        Return the matrix representations from the database.

        INPUT:

        - ``representation_type`` -- instance of enum
          :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation


        OUTPUT:

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: CHA3 = algebras.CubicHecke(2)
            sage: GER = CHA3.extension_ring(generic=True)
            sage: cha_db = CHA3._database
            sage: rt = CHA3.repr_type
            sage: m1 =cha_db.read_matrix_representation(rt.SplitIrredMarin, 1, 3, GER)
            sage: len(m1)
            7
            sage: GBR = CHA3.base_ring(generic=True)
            sage: m1rl = cha_db.read_matrix_representation(rt.RegularLeft, 1, 3, GBR)
            sage: m1rl[0].dimensions()
            (24, 24)
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType, GenSign
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' %(RepresentationType))

        td = ring_of_definition.gens_dict_recursive()
        if 'e3' in td.keys():
            td['j'] = td['e3']

        num_rep = representation_type.number_of_representations(nstrands)
        rep_list = self.read(representation_type.data_section(), translation_dict=td, nstrands=nstrands)
        if gen_ind > 0 :
            rep_list = [rep_list[GenSign.pos][i] for i in range(num_rep)]
            matrix_list = [matrix(ring_of_definition, rep[gen_ind-1], sparse=True) for rep in rep_list]
        else:
            # data of inverse of generators is stored under negative strand-index
            rep_list = [rep_list[GenSign.neg][i] for i in range(num_rep) ]
            matrix_list = [matrix(ring_of_definition, rep[-gen_ind-1], sparse=True) for rep in rep_list]
        for m in matrix_list: m.set_immutable()
        return matrix_list


class CubicHeckeFileCache(SageObject):
    """
    A class to cache calculations of the :class:`CubicHeckeAlgebra` in the local
    file system.
    """

    class section(Enum):
        r"""
        Enum for the different sections of file cache. The following choices are
        possible:

        - ``matrix_representations``  -- file cache for representation matrices
          of basis elements
        - ``braid_images``  -- file cache for images of braids
        - ``basis_extensions`` -- file cache for a dynamical growing basis used
          in the case of cubic Hecke algebras on more than 4 strands
        - ``markov_trace`` -- file cache for intermediate results of long
          calculations in order to recover the results already obtained by
          preboius attemps of calculation until the corresponding intermediate
          step

        Examples::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(CHA2)
            sage: cha_fc.section
            <enum 'section'>
        """

        def filename(self, nstrands=None):
            r"""
            Return the file name under which the data of this file cache section
            is stored as an sobj-file.

            INPUT:

            - ``nstrands`` -- Integer, number of strands of the underlying braid
              group if the data file depends on it. Otherwise use default ``None``

            Examples::

                sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
                sage: CHA2 = algebras.CubicHecke(2)
                sage: cha_fc = CubicHeckeFileCache(CHA2)
                sage: cha_fc.section.matrix_representations.filename(2)
                'matrix_representations_2.sobj'
                sage: cha_fc.section.braid_images.filename(2)
                'braid_images_2.sobj'
            """
            if nstrands is None:
                return '%s.sobj' %(self.value)
            else:
                return '%s_%s.sobj' %(self.value, nstrands)

        matrix_representations  = 'matrix_representations'
        braid_images            = 'braid_images'
        basis_extensions        = 'basis_extensions'
        markov_trace            = 'markov_trace'


    def __init__(self, num_strands):
        r"""
        Python constructor.

        INPUT:

        - ``cubic_hecke_algebra`` -- instance of :class:`CubicHeckeAlgebra`
          whose data should be cached in the file system.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc._file_cache_path.endswith('cubic_hecke')
            True
        """
        self._nstrands      = num_strands

        from sage.env import DOT_SAGE
        self._file_cache_path = os.path.join(DOT_SAGE, 'cubic_hecke')
        self._data_library = {}

        from sage.misc.misc import sage_makedirs
        sage_makedirs(self._file_cache_path)

    def reset_library(self, section=None):
        r"""
        Reset the file cache corresponding to the specified ``section``.

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
            sage: cha2_fc.reset_library(cha2_fc.section.matrix_representations)
            sage: data_mat = cha2_fc.read(cha2_fc.section.matrix_representations)
            sage: len(data_mat.keys())
            4
        """
        if section is None:
            for sec in self.section:
                self.reset_library(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' %(CubicHeckeFileCache.section))
  
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        data_lib = self._data_library
        empty_dict = {}
        if  section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                new_dict={}
                empty_dict.update({rep_type:new_dict})
        elif  section == self.section.basis_extensions:
            empty_dict = []
        data_lib.update({section:empty_dict})


    def is_empty(self, section=None):
        r"""
        Return ``True`` if the cache of the given ``section`` is empty. 

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections
          

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library()
            sage: cha2_fc.is_empty()
            True
        """
        if section is None:
            return all(self.is_empty(section=sec) for sec in self.section)

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' %(CubicHeckeFileCache.section))
  
        self.read(section)
        data_lib = self._data_library[section]
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if  section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                if len(data_lib[rep_type]) > 0:
                    return False
            return True

        if  section == self.section.basis_extensions and self._nstrands > 4:
            # the new generators and their inverses are not counted
            # since they are added during initialization
            return len(data_lib) <= 2*(self._nstrands -4)
        return len(data_lib) == 0
       


    # --------------------------------------------------------------------------
    # save data file system
    # --------------------------------------------------------------------------
    def write(self, section=None):
        r"""
        Write data from memory to the file system.

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to.
          If omitted data of all sections is written to the file system
        - ``step`` -- integer, to indicate the intermediate step in the
          calculation of the Markov trace coefficients. This makes sence for
          the ``markov_trace`` section, only

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.write(cha2_fc.section.braid_images)
        """
        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section is None:
            for sec in self.section:
                if sec in data_lib.keys():
                   self.write(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' %(CubicHeckeFileCache.section))

        if section not in data_lib.keys():
            raise ValueError("No data for file %s in memory" %(section))

        verbose('saving file cache %s ...' %(section))
        fname = os.path.join(lib_path, section.filename(self._nstrands))
        with atomic_write(fname, binary=True) as f:
            f.write(_base_dumps(data_lib[section]))
            f.close()

    # --------------------------------------------------------------------------
    # read from file system
    # --------------------------------------------------------------------------
    def read(self, section):
        r"""
        Read data into memory from the file system.
        

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to

        OUTPUT:

        Dictionary containing the data library corresponding to the section
        of file cache
        
        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
        """
        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' %(CubicHeckeFileCache.section))

        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section in data_lib.keys():
            return data_lib[section]

        verbose('loading file cache %s ...' %(section))
        fname = os.path.join(lib_path, section.filename(self._nstrands))
        try:
            data_lib[section] = load(fname)
            verbose('... finished!')
        except IOError:
            self.reset_library(section)
            verbose('... not found!')

        return data_lib[section]




    # --------------------------------------------------------------------------
    # matrix_reprs_from_file_cache_
    # --------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, monomial_tietze, ring_of_definition):
        r"""
        Return the matrix representations of the given monomial (in Tietze form)
        if it has been stored in the file cache before. Otherwise ``None`` is
        returned.

        INPUT:

        - ``representation_type`` -- instance of enum :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation

        - ``monomial_tietze`` -- tuple representing the braid in Tietze form

        - ``ring_of_definition`` -- instance of :class:`CubicHeckeRingOfDefinition` resp.
          :class:`CubicHeckeExtensionRing` (depending whether ``representation_type``
          is split or not)

        OUTPUT:

        Dictionary containing all matrix representations of ``self`` of the given
        ``representation_type`` which have been stored in the file cache.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gt = g.Tietze()
            sage: rt = CHA2.repr_type
            sage: g.matrix(representation_type=rt.RegularLeft)
            [ 0 -v  1]
            [ 1  u  0]
            [ 0  w  0]
            sage: [_] == cha_fc.read_matrix_representation(rt.RegularLeft, gt, R)
            True
            sage: cha_fc.reset_library(cha_fc.section.matrix_representations)
            sage: cha_fc.write(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, gt, R) == None
            True
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' %(RepresentationType))

        matrix_representations = self.read(self.section.matrix_representations)[representation_type]
        if monomial_tietze in matrix_representations.keys():
            matrix_list_dict = matrix_representations[monomial_tietze]
            matrix_list = [matrix(ring_of_definition, mat_dict, sparse=True) for mat_dict in matrix_list_dict]
            for m in matrix_list: m.set_immutable()
            return matrix_list
        return None




    # --------------------------------------------------------------------------
    # matrix_representation to file cache
    # --------------------------------------------------------------------------
    def write_matrix_representation(self, representation_type, monomial_tietze, matrix_list):
        r"""
        Write the matrix representation of a monomial to the file cache.

        INPUT:

        - ``representation_type`` -- instance of enum :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation

        - ``monomial_tietze`` -- tuple representing the braid in Tietze form

        - ``matrix_list`` -- list of matrices corresponding to the irreducible
          representations

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gi = ~g; git = gi.Tietze()
            sage: rt = CHA2.repr_type
            sage: m = gi.matrix(representation_type=rt.RegularRight)
            sage: cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            [
            [     0      1 (-u)/w]
            [     0      0    1/w]
            [     1      0    v/w]
            ]
            sage: CHA2.reset_filecache(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, git, R) == None
            True
            sage: cha_fc.write_matrix_representation(rt.RegularRight, git, [m])
            sage: [m] == cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            True
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' %(RepresentationType))

        matrix_representations = self.read(self.section.matrix_representations)[representation_type]

        if monomial_tietze in matrix_representations.keys():
            # entry already registered
            return

        matrix_representation_dict = [simplify(mat) for mat in list(matrix_list)]
        matrix_representations[monomial_tietze] = matrix_representation_dict

        self.write(self.section.matrix_representations)
        return

    # --------------------------------------------------------------------------
    # read braid images from file cache
    # --------------------------------------------------------------------------
    def read_braid_image(self, braid_tietze, ring_of_definition):
        r"""
        Return the list of pre calculated braid images from file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form

        - ``ring_of_definition`` -- instance of :class:`CubicHeckeRingOfDefinition`

        OUTPUT:

        A dictionary containing the pre calculated braid image of the given
        braid.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b2 = b**2
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
            sage: b2_img = CHA2(b2); b2_img
            (-v) + u*c + w*c^-1
            sage: cha_fc.write_braid_image(b2.Tietze(), b2_img.to_vector())
            sage: cha_fc.read_braid_image(b2.Tietze(), ring_of_definition)
            (-v, u, w)
        """
        braid_images = self.read(self.section.braid_images)
        if braid_tietze in braid_images.keys():
            braid_image = braid_images[braid_tietze]
            result_list = [ring_of_definition(cf) for cf in list(braid_image)]
            from sage.modules.free_module_element import vector
            return vector(ring_of_definition, result_list)
        return None





    # --------------------------------------------------------------------------
    # braid image to_file cache
    # --------------------------------------------------------------------------
    def write_braid_image(self, braid_tietze, braid_image_vect):
        r"""
        Write the braid image of the given braid to the file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form
        - ``braid_image_vect`` -- image of the given braid as a vector with
          respect to the basis of the cubic Hecke algebra

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b3 = b**3
            sage: b3_img = CHA2(b3); b3_img
            (-u*v+w) + (u^2-v)*c + u*w*c^-1
            sage: cha_fc.write_braid_image(b3.Tietze(), b3_img.to_vector())
            sage: cha_fc.read_braid_image(b3.Tietze(), ring_of_definition)
            (-u*v + w, u^2 - v, u*w)
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.write(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
        """
        braid_images = self.read(self.section.braid_images)

        if braid_tietze in braid_images.keys():
            # entry already registered
            return

        braid_image_dict = [simplify(cf) for cf in list(braid_image_vect)]
        braid_images[braid_tietze] = braid_image_dict

        self.write(self.section.braid_images)
        return


    # --------------------------------------------------------------------------
    # basis to file cache
    # --------------------------------------------------------------------------
    def update_basis_extensions(self, new_basis_extensions):
        r"""
        Update the file cache for basis extensions for cubic Hecke algebras on
        more than 4 strands according to the given ``new_basis_extensions``.

        INPUT:

        - ``new_basis_extensions`` -- list of additional (to the static basis)
          basis elements which should replace the former such list in the file.
 
        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
            sage: cha_fc.update_basis_extensions([(1,), (-1,)])
            sage: cha_fc.read(CubicHeckeFileCache.section.basis_extensions)
            [(1,), (-1,)]
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.write(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
        """
        self._data_library.update({self.section.basis_extensions:new_basis_extensions})
        self.write(self.section.basis_extensions)
        return

    # --------------------------------------------------------------------------
    # Intermediate results of Markov trace coefficient calculation
    # --------------------------------------------------------------------------
    def markov_trace(self, step, target=None):
        r"""
        Return a wrapper to store intermediate results obtained during the calculation
        of the Markov trace coefficients.

        INPUT:

        - ``step`` -- integer or string, to indicate the intermediate step in the
          calculation of the Markov trace coefficients. This makes sence for
          the ``markov_trace`` section, only
        """

        class wrapper:
            r"""
            Wrapper that stores intermediate results obtained during the calculation
            of the Markov trace coefficients.

            INPUT (to the constructor):

            - ``filecache`` -- pointer to the outer class
            - ``step`` -- according to the method
            """
            def __init__(self, filecache, step, target=None):
                r"""
                """
                self._fc = filecache
                self._step = step
                self._target = target
                self._target_key = 'targets'
                self._data_section = None
                self._data_step = None
                self._message = 'step %s for Markov trace coefficients on %s strands' %(step, filecache._nstrands)
                self._time = verbose('Starting calulation of %s' %(self._message))

            def _set_data_step(self, data):
                r"""
                """
                if type(data) is list:
                    from sage.rings.localization import LocalizationElement
                    if isinstance(data[0], LocalizationElement):
                         data = [simplify(item) for item in data]
                self._data_step = data
                self._time = verbose('Calculation finished for %s' %(self._message), t=self._time)
                return data

            def _set_target(self):
                r"""
                Store the target of this step in the section space.
                """
                if not self._target:
                    return

                if not self._data_section:
                    return

                step = self._step
                target = self._target
                target_key = self._target_key
                if target_key in self._data_section:
                    targets = self._data_section[target_key]
                else:
                    targets = {}
                    self._data_section[target_key] = targets

                if  target in targets.keys():
                    step_list = targets[target]
                    if not step in step_list:
                        step_list.append(step)
                else:
                    targets[target] = [step]
                self._time = verbose('Target %s set for %s' %(target, self._message), t=self._time)

            def remove_temporary_data(self):
                r"""
                Remove data of temporary steps if target is complete.
                """
                fc   = self._fc
                step = self._step
                target_key = self._target_key
                data_section = self._data_section

                if not data_section:
                    return

                if target_key not in data_section:
                    return

                targets = data_section[target_key]
                if step not in targets.keys():
                    return

                for sub_step in targets[step]:
                    self._time = verbose('Remove temporary step %s set for %s' %(sub_step, self._message), t=self._time)
                    data_section.pop(sub_step)

                fc.write(fc.section.markov_trace)

            def __enter__(self):
                r"""
                OUTPUT:

                    A pair (data, method) one of which is ``None``.
                """
                fc   = self._fc
                step = self._step
                data_section = fc.read(fc.section.markov_trace)
                self._data_section = data_section
                if self._target:
                    self._set_target()
                if step in data_section.keys():
                    self._data_step = data_section[step]
                    self._time = verbose('Found previous result for %s' %(self._message), t=self._time)
                    return self._data_step, self._set_data_step
                return None, self._set_data_step

            def __exit__(self, exc_type, exc_val, exc_tb):
                r"""
                """
                fc   = self._fc
                step = self._step
                self._data_section[step] = self._data_step
                fc.write(fc.section.markov_trace)

        return wrapper(self, step, target)

# -----------------------------------------------------------------------------
# Demo data section
# -----------------------------------------------------------------------------
template=r"""# generated demo data for cubic Hecke algebra
def read_%s_generated(%snstands):
    r\"""
    %s
    \"""
    data = {}
    data[2] = %s
    data[3] = %s
    return data[nstands]

"""

doc=r"""Return precomputed data of Ivan Marin

    INPUT:

    - ``a, b, c, j, u, v, w`` -- variables for the representation matrices
    - ``nstarnds`` -- integer, number of strands of the cubic Hecke algebra in
      question

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_bas_generated
        sage: read_bas_generated(2)
        [[], [1], [-1]]
"""


def create_demo_data(filename='demo_data.py'):
        r"""
        Return the input for the statical ``data_library`` as a dictionary
        such that it can be stored in this file.

        This makes only sense if the feature is present.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import create_demo_data
            sage: create_demo_data()
        """
        from sage.algebras.hecke_algebras.cubic_hecke_algebra import CubicHeckeAlgebra
        from textwrap import fill
        def fl(s):
            return fill(str(s), subsequent_indent='        ')
        from sympy import var
        TE = {str(x):x for x in var('a, b, c, j')}
        TB = {str(x):x for x in var('u, v, w')}

        from database_cubic_hecke import read_basis, read_irr, read_reg
        bas2  = fl(read_basis(num_strands=2))
        bas3  = fl(read_basis(num_strands=3))
        irr2  = fl(read_irr(translation_dict=TE, num_strands=2))
        irr3  = fl(read_irr(translation_dict=TE, num_strands=3))
        regl2 = fl(read_reg(translation_dict=TB, num_strands=2))
        regl3 = fl(read_reg(translation_dict=TB, num_strands=3))
        regr2 = fl(read_reg(translation_dict=TB, num_strands=2, right=True))
        regr3 = fl(read_reg(translation_dict=TB, num_strands=3, right=True))

        bas  = template %('bas', '', doc, bas2, bas3)
        irr  = template %('irr', 'a, b, c, j, ', doc, irr2, irr3)
        regl = template %('regl', 'u, v, w, ', doc, regl2, regl3)
        regr = template %('regr', 'u, v, w, ', doc, regr2, regr3)
        with open(filename, 'w') as f:
            f.write(bas)
            f.write(irr)
            f.write(regl)
            f.write(regr)


# generated demo data for cubic Hecke algebra
def read_bas_generated(nstands):
    r"""
    Return precomputed data of Ivan Marin

    INPUT:

    - ``a, b, c, j, u, v, w`` -- variables for the representation matrices
    - ``nstarnds`` -- integer, number of strands of the cubic Hecke algebra in
      question

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_bas_generated
        sage: read_bas_generated(2)
        [[], [1], [-1]]

    """
    data = {}
    data[2] = [[], [1], [-1]]
    data[3] = [[], [1], [-1], [2], [-2], [1, 2], [1, -2], [-1, 2], [-1, -2], [1, 2,
        1], [1, 2, -1], [-1, 2, 1], [-1, 2, -1], [1, -2, 1], [-1, -2,
        1], [2, 1], [-2, 1], [2, -1], [-2, -1], [1, -2, -1], [-1, -2,
        -1], [2, -1, 2], [1, 2, -1, 2], [-1, 2, -1, 2]]
    return data[nstands]

# generated demo data for cubic Hecke algebra
def read_irr_generated(a, b, c, j, nstands):
    r"""
    Return precomputed data of Ivan Marin

    INPUT:

    - ``a, b, c, j, u, v, w`` -- variables for the representation matrices
    - ``nstarnds`` -- integer, number of strands of the cubic Hecke algebra in
      question

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_irr_generated
        sage: read_irr_generated(3, 5, 7, 9, 2)
        ([1, 1, 1],
         [[{(0, 0): 3}], [{(0, 0): 7}], [{(0, 0): 5}]],
         [[{(0, 0): 1/3}], [{(0, 0): 1/7}], [{(0, 0): 1/5}]])

    """
    data = {}
    data[2] = ([1, 1, 1], [[{(0, 0): a}], [{(0, 0): c}], [{(0, 0): b}]], [[{(0, 0):
        1/a}], [{(0, 0): 1/c}], [{(0, 0): 1/b}]])
    data[3] = ([1, 1, 1, 2, 2, 2, 3], [[{(0, 0): a}, {(0, 0): a}], [{(0, 0): c},
        {(0, 0): c}], [{(0, 0): b}, {(0, 0): b}], [{(0, 0): b, (1, 0):
        b*c, (1, 1): c}, {(0, 0): c, (0, 1): -1, (1, 1): b}], [{(0,
        0): a, (1, 0): a*b, (1, 1): b}, {(0, 0): b, (0, 1): -1, (1,
        1): a}], [{(0, 0): a, (1, 0): a*c, (1, 1): c}, {(0, 0): c, (0,
        1): -1, (1, 1): a}], [{(0, 0): c, (1, 0): a*c + b**2, (1, 1):
        b, (2, 0): b, (2, 1): 1, (2, 2): a}, {(0, 0): a, (0, 1): -1,
        (0, 2): b, (1, 1): b, (1, 2): -a*c - b**2, (2, 2): c}]],
        [[{(0, 0): 1/a}, {(0, 0): 1/a}], [{(0, 0): 1/c}, {(0, 0):
        1/c}], [{(0, 0): 1/b}, {(0, 0): 1/b}], [{(0, 0): 1/b, (1, 0):
        -1, (1, 1): 1/c}, {(0, 0): 1/c, (0, 1): 1/(b*c), (1, 1):
        1/b}], [{(0, 0): 1/a, (1, 0): -1, (1, 1): 1/b}, {(0, 0): 1/b,
        (0, 1): 1/(a*b), (1, 1): 1/a}], [{(0, 0): 1/a, (1, 0): -1, (1,
        1): 1/c}, {(0, 0): 1/c, (0, 1): 1/(a*c), (1, 1): 1/a}], [{(0,
        0): 1/c, (1, 0): -a/b - b/c, (1, 1): 1/b, (2, 0): 1/b, (2, 1):
        -1/(a*b), (2, 2): 1/a}, {(0, 0): 1/a, (0, 1): 1/(a*b), (0, 2):
        1/b, (1, 1): 1/b, (1, 2): a/b + b/c, (2, 2): 1/c}]])
    return data[nstands]

# generated demo data for cubic Hecke algebra
def read_regl_generated(u, v, w, nstands):
    r"""
    Return precomputed data of Ivan Marin

    INPUT:

    - ``a, b, c, j, u, v, w`` -- variables for the representation matrices
    - ``nstarnds`` -- integer, number of strands of the cubic Hecke algebra in
      question

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_regl_generated
        sage: read_regl_generated(3, 5, 7, 2)
        ([3],
         [[{(0, 1): -5, (0, 2): 1, (1, 0): 1, (1, 1): 3, (2, 1): 7}]],
         [[{(0, 1): 1, (0, 2): -3/7, (1, 2): 1/7, (2, 0): 1, (2, 2): 5/7}]])

    """
    data = {}
    data[2] = ([3], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2):
        v/w}]])
    data[3] = ([24], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w, (3,
        5): -v, (3, 7): 1, (4, 6): -v, (4, 8): 1, (5, 3): 1, (5, 5):
        u, (6, 4): 1, (6, 6): u, (7, 5): w, (8, 6): w, (9, 9): u, (9,
        15): 1, (10, 10): u, (10, 17): 1, (11, 9): w, (12, 10): w,
        (13, 13): u, (13, 16): 1, (14, 13): w, (15, 9): -v, (15, 11):
        1, (16, 13): -v, (16, 14): 1, (17, 10): -v, (17, 12): 1, (18,
        19): -v, (18, 20): 1, (19, 18): 1, (19, 19): u, (20, 19): w,
        (21, 22): -v, (21, 23): 1, (22, 21): 1, (22, 22): u, (23, 22):
        w}, {(0, 3): -v, (0, 4): 1, (1, 15): -v, (1, 16): 1, (1, 22):
        -v, (1, 23): u*v/w, (2, 17): -v, (2, 18): 1, (2, 23): v*(u*v -
        w)/w, (3, 0): 1, (3, 3): u, (4, 3): w, (5, 9): -v, (5, 10): 1,
        (5, 12): -u/w, (5, 22): u, (5, 23): -u**2/w, (6, 11): -v, (6,
        22): w, (6, 23): -u, (7, 12): -u*v/w, (7, 13): -v, (7, 19): 1,
        (7, 21): -v, (7, 23): -u**2*v/w, (8, 12): -v, (8, 14): -v, (8,
        20): 1, (8, 23): -u*v, (9, 5): 1, (9, 9): u, (9, 23): u/w,
        (10, 9): w, (10, 11): -u, (11, 6): 1, (11, 11): u, (11, 13):
        u, (12, 13): w, (12, 23): -v, (13, 23): 1, (14, 8): 1, (14,
        14): u, (14, 23): v, (15, 1): 1, (15, 12): u/w, (15, 15): u,
        (16, 11): v, (16, 15): w, (16, 23): -u, (17, 2): 1, (17, 12):
        u*v/w, (17, 17): u, (18, 12): v, (18, 17): w, (19, 21): w,
        (19, 23): v, (20, 14): w, (21, 7): 1, (21, 21): u, (21, 23):
        u*v/w, (22, 11): 1, (23, 12): 1, (23, 23): u}]], [[{(0, 1): 1,
        (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2): v/w, (3, 5): 1,
        (3, 7): -u/w, (4, 6): 1, (4, 8): -u/w, (5, 7): 1/w, (6, 8):
        1/w, (7, 3): 1, (7, 7): v/w, (8, 4): 1, (8, 8): v/w, (9, 11):
        1/w, (10, 12): 1/w, (11, 11): v/w, (11, 15): 1, (12, 12): v/w,
        (12, 17): 1, (13, 14): 1/w, (14, 14): v/w, (14, 16): 1, (15,
        9): 1, (15, 11): -u/w, (16, 13): 1, (16, 14): -u/w, (17, 10):
        1, (17, 12): -u/w, (18, 19): 1, (18, 20): -u/w, (19, 20): 1/w,
        (20, 18): 1, (20, 20): v/w, (21, 22): 1, (21, 23): -u/w, (22,
        23): 1/w, (23, 21): 1, (23, 23): v/w}, {(0, 3): 1, (0, 4):
        -u/w, (1, 15): 1, (1, 16): -u/w, (1, 22): u*v/w, (1, 23):
        -u/w, (2, 17): 1, (2, 18): -u/w, (3, 4): 1/w, (4, 0): 1, (4,
        4): v/w, (5, 9): 1, (5, 10): -u/w, (5, 13): -u/w, (5, 22):
        -u**2/w, (6, 11): 1, (6, 12): -u/w, (6, 13): -u*v/w, (6, 22):
        -u, (7, 19): -u/w, (7, 21): 1, (8, 13): -v, (8, 14): 1, (8,
        20): -u/w, (9, 10): 1/w, (9, 22): u/w, (10, 5): 1, (10, 6):
        -u/w, (10, 10): v/w, (10, 13): -u**2/w, (10, 23): u/w, (11,
        22): 1, (12, 13): -u, (12, 23): 1, (13, 12): 1/w, (13, 13):
        v/w, (14, 20): 1/w, (15, 13): u/w, (15, 16): 1/w, (15, 22):
        -v/w, (16, 1): 1, (16, 6): v/w, (16, 13): u*v/w, (16, 16):
        v/w, (17, 13): u*v/w, (17, 18): 1/w, (17, 23): -v/w, (18, 2):
        1, (18, 13): v, (18, 18): v/w, (18, 23): -v**2/w, (19, 7): 1,
        (19, 12): v/w, (19, 19): v/w, (19, 23): u*v/w, (20, 8): 1,
        (20, 20): v/w, (20, 23): v, (21, 13): -v/w, (21, 19): 1/w,
        (22, 6): 1/w, (22, 13): u/w, (22, 22): v/w, (23, 13): 1}]])
    return data[nstands]

# generated demo data for cubic Hecke algebra
def read_regr_generated(u, v, w, nstands):
    r"""
    Return precomputed data of Ivan Marin

    INPUT:

    - ``a, b, c, j, u, v, w`` -- variables for the representation matrices
    - ``nstarnds`` -- integer, number of strands of the cubic Hecke algebra in
      question

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import read_regr_generated
        sage: read_regr_generated(3, 5, 7, 2)
        ([3],
         [[{(0, 1): -5, (0, 2): 1, (1, 0): 1, (1, 1): 3, (2, 1): 7}]],
         [[{(0, 1): 1, (0, 2): -3/7, (1, 2): 1/7, (2, 0): 1, (2, 2): 5/7}]])

    """
    data = {}
    data[2] = ([3], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w}]],
        [[{(0, 1): 1, (0, 2): -u/w, (1, 2): 1/w, (2, 0): 1, (2, 2):
        v/w}]])
    data[3] = ([24], [[{(0, 1): -v, (0, 2): 1, (1, 0): 1, (1, 1): u, (2, 1): w, (3,
        15): -v, (3, 17): 1, (4, 16): -v, (4, 18): 1, (4, 22): v**2,
        (4, 23): -v, (5, 9): -v, (5, 10): 1, (6, 13): -v, (6, 19): 1,
        (6, 21): -v, (6, 22): -u*v, (7, 11): -v, (7, 12): 1, (8, 14):
        -v, (8, 20): 1, (8, 22): -v*w, (9, 5): 1, (9, 9): u, (9, 23):
        u/w, (10, 9): w, (10, 21): -u, (10, 22): -u**2, (11, 7): 1,
        (11, 11): u, (11, 21): u, (11, 23): u*v/w, (12, 11): w, (12,
        22): -u*w, (13, 6): 1, (13, 13): u, (13, 22): v, (14, 8): 1,
        (14, 14): u, (14, 23): v, (15, 3): 1, (15, 15): u, (15, 22):
        u, (15, 23): -u**2/w, (16, 4): 1, (16, 16): u, (16, 21): v,
        (17, 15): w, (17, 22): u*v, (17, 23): -u, (18, 16): w, (19,
        13): w, (20, 14): w, (21, 22): -v, (21, 23): 1, (22, 21): 1,
        (22, 22): u, (23, 22): w}, {(0, 3): -v, (0, 4): 1, (1, 5): -v,
        (1, 6): 1, (2, 7): -v, (2, 8): 1, (3, 0): 1, (3, 3): u, (4,
        3): w, (5, 1): 1, (5, 5): u, (6, 5): w, (7, 2): 1, (7, 7): u,
        (8, 7): w, (9, 9): u, (9, 15): 1, (10, 13): u, (10, 16): 1,
        (10, 22): -v, (11, 9): w, (12, 13): w, (12, 23): -v, (13, 23):
        1, (14, 21): w, (14, 23): v, (15, 9): -v, (15, 11): 1, (16,
        22): w, (16, 23): -u, (17, 13): -v, (17, 14): 1, (17, 21): -v,
        (18, 19): -v, (18, 20): 1, (19, 18): 1, (19, 19): u, (20, 19):
        w, (21, 17): 1, (21, 21): u, (22, 10): 1, (22, 22): u, (23,
        12): 1, (23, 23): u}]], [[{(0, 1): 1, (0, 2): -u/w, (1, 2):
        1/w, (2, 0): 1, (2, 2): v/w, (3, 15): 1, (3, 17): -u/w, (3,
        23): u*(u*v - w)/w**2, (4, 16): 1, (4, 18): -u/w, (4, 22): -v,
        (4, 23): u*v/w, (5, 9): 1, (5, 10): -u/w, (5, 21): -u/w, (5,
        22): -u**2/w, (5, 23): -u*v/w**2, (6, 13): 1, (6, 19): -u/w,
        (6, 23): -v/w, (7, 11): 1, (7, 12): -u/w, (7, 21): -u*v/w, (7,
        22): -u, (7, 23): -u*v**2/w**2, (8, 14): 1, (8, 20): -u/w, (8,
        21): -v, (8, 23): -v**2/w, (9, 10): 1/w, (9, 22): u/w, (10,
        5): 1, (10, 10): v/w, (10, 22): u*v/w, (11, 12): 1/w, (11,
        23): u/w, (12, 7): 1, (12, 12): v/w, (12, 23): u*v/w, (13,
        19): 1/w, (14, 20): 1/w, (15, 17): 1/w, (15, 21): u/w, (16,
        18): 1/w, (17, 3): 1, (17, 17): v/w, (17, 21): u*v/w, (18, 4):
        1, (18, 18): v/w, (18, 21): v, (19, 6): 1, (19, 19): v/w, (19,
        22): v, (20, 8): 1, (20, 20): v/w, (20, 23): v, (21, 22): 1,
        (21, 23): -u/w, (22, 23): 1/w, (23, 21): 1, (23, 23): v/w},
        {(0, 3): 1, (0, 4): -u/w, (1, 5): 1, (1, 6): -u/w, (2, 7): 1,
        (2, 8): -u/w, (3, 4): 1/w, (4, 0): 1, (4, 4): v/w, (5, 6):
        1/w, (6, 1): 1, (6, 6): v/w, (7, 8): 1/w, (8, 2): 1, (8, 8):
        v/w, (9, 11): 1/w, (10, 13): -u**2/w, (10, 16): -u/w, (10,
        22): 1, (11, 11): v/w, (11, 15): 1, (12, 13): -u, (12, 23): 1,
        (13, 12): 1/w, (13, 13): v/w, (14, 12): v/w, (14, 14): v/w,
        (14, 17): 1, (15, 9): 1, (15, 11): -u/w, (16, 10): 1, (16,
        12): -u/w, (16, 16): v/w, (17, 13): u*v/w, (17, 14): -u/w,
        (17, 21): 1, (18, 19): 1, (18, 20): -u/w, (19, 20): 1/w, (20,
        18): 1, (20, 20): v/w, (21, 13): -v/w, (21, 14): 1/w, (22,
        13): u/w, (22, 16): 1/w, (23, 13): 1}]])
    return data[nstands]

