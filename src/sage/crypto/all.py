from __future__ import absolute_import
import sage.crypto.sbox
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.crypto.mq.sbox', 'SBox', sage.crypto.sbox.SBox)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.classical', ['AffineCryptosystem',
                                      'HillCryptosystem',
                                      'SubstitutionCryptosystem',
                                      'ShiftCryptosystem',
                                      'TranspositionCryptosystem',
                                      'VigenereCryptosystem',
                                     ])

lazy_import('sage.crypto.stream', ['LFSRCryptosystem',
                                   'ShrinkingGeneratorCryptosystem',
                                  ])

lazy_import('sage.crypto.lfsr', ['lfsr_sequence',
                                 'lfsr_autocorrelation',
                                 'lfsr_connection_polynomial',
                                ])
