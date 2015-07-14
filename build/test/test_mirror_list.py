# -*- coding: utf-8 -*-
"""
Test downloading files
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import unittest
import tempfile

from .capture import CapturedLog
from sage_bootstrap.mirror_list import MirrorList



class MirrorListTestCase(unittest.TestCase):

    def test_mirror_list(self):
        with CapturedLog() as log:
            ml = MirrorList()
        msg = log.messages()
        if len(msg) > 0:
            self.assertEqual(msg[0], ('INFO', 'Downloading the Sage mirror list'))
        self.assertTrue(len(ml.mirrors) >= 0)
        self.assertTrue(ml.fastest.startswith('http://'))
