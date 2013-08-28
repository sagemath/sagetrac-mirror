# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
cryptodir = $(pythondir)/sage/crypto
crypto_PYTHON = lattice.py __init__.py stream_cipher.py cryptosystem.py stream.py all.py classical.py lwe.py classical_cipher.py cipher.py lfsr.py util.py

AM_DEFAULT_SOURCE_EXT = .pyx
crypto_LTLIBRARIES = boolean_function.la

crypto_mqdir = $(cryptodir)/mq
crypto_mq_PYTHON = \
    mq/__init__.py \
    mq/mpolynomialsystemgenerator.py \
    mq/sbox.py \
    mq/sr.py

block_cipherdir = $(cryptodir)/block_cipher
block_cipher_PYTHON = \
    block_cipher/all.py \
    block_cipher/__init__.py \
    block_cipher/miniaes.py \
    block_cipher/sdes.py

crypto_public_keydir = $(cryptodir)/public_key
crypto_public_key_PYTHON = \
    public_key/all.py \
    public_key/blum_goldwasser.py \
    public_key/__init__.py
