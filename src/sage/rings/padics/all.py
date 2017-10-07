from __future__ import absolute_import
from .generic_nodes import is_pAdicField, is_pAdicRing
from .factory import Zp, Zq, Zp as pAdicRing, ZpCR, ZpCA, ZpFM, ZpFP, ZqCR, ZqCA, ZqFM, ZqFP #, ZpL, ZqL
from .factory import Qp, Qq, Qp as pAdicField, QpCR, QpFP, QqCR, QqFP #, QpL, QqL
from .factory import pAdicExtension
from .padic_generic import local_print_mode
from .pow_computer import PowComputer
from .pow_computer_ext import PowComputer_ext_maker
from .discrete_value_group import DiscreteValueGroup
