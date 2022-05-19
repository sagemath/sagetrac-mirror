# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``imagemagick``

Currently we only check for the presence of ``convert``. When needed other
commands like ``magick``, ``magick-script``, ``convert``, ``mogrify``,
``identify``, ``composite``, ``montage``, ``compare``, etc. could be also
checked in this module.
"""
# ****************************************************************************
#       Copyright (C) 2018 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable, FeatureTestResult
from .join_feature import JoinFeature

class Convert(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``convert``

    EXAMPLES::

        sage: from sage.features.imagemagick import Convert
        sage: Convert().is_present()  # optional - imagemagick
        FeatureTestResult('convert', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.imagemagick import Convert
            sage: isinstance(Convert(), Convert)
            True
        """
        Executable.__init__(self, "convert", executable="convert")

    def is_functional(self):
        r"""
        Return whether command ``convert`` in the path is functional.

        EXAMPLES::

            sage: from sage.features.imagemagick import Convert
            sage: Convert().is_functional()   # optional - imagemagick
            FeatureTestResult('convert', True)

        """
        # Create the content of 1-pixel png file
        content = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x00\x00\x00\x00:~\x9bU\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01H\xaf\xa4q\x00\x00\x00\x00IEND\xaeB`\x82'

        # NOTE:
        #
        # This is how the above content of a 1 pixel PNG was created::
        #
        #    sage: import numpy as np
        #    sage: from PIL import Image
        #    sage: image = Image.fromarray(np.array([[100]], dtype=np.uint8))
        #    sage: image.save('file.png')
        #    sage: with open('file.png', 'rb') as f:
        #    ....:     content = f.read()

        # create a png file with the content
        from sage.misc.temporary_file import tmp_filename
        base_filename_png = tmp_filename(ext='.png')
        with open(base_filename_png, 'wb') as f:
            f.write(content)

        # Set up filenames
        import os
        base, filename_png = os.path.split(base_filename_png)
        filename, _png = os.path.splitext(filename_png)
        filename_gif = filename + '.gif'

        # running command convert (taken from sage/plot/animate.py)
        from subprocess import run
        cmd = ['convert', '-dispose', 'Background', '-delay', '20',
                '-loop', '0', filename_png, filename_gif]

        try:
            result = run(cmd, cwd=base, capture_output=True, text=True)
        except OSError as e:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'raised an OSError "{}" '.format(' '.join(cmd), e))

        # If an error occurred, return False
        if result.returncode:
            return FeatureTestResult(self, False, reason='Running command "{}" '
                        'returned non-zero exit status "{}" with stderr '
                        '"{}" and stdout "{}".'.format(result.args,
                                                       result.returncode,
                                                       result.stderr.strip(),
                                                       result.stdout.strip()))

        # If necessary, run more tests here
        # ...

        # The command seems functional
        return FeatureTestResult(self, True)

class ImageMagick(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    ``ImageMagick``

    Currently, only the availability of ``convert`` is checked.

    EXAMPLES::

        sage: from sage.features.imagemagick import ImageMagick
        sage: ImageMagick().is_present()  # optional - imagemagick
        FeatureTestResult('imagemagick', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.imagemagick import ImageMagick
            sage: isinstance(ImageMagick(), ImageMagick)
            True
        """
        JoinFeature.__init__(self, "imagemagick",
                             [Convert()],
                             spkg="imagemagick",
                             url="https://www.imagemagick.org/")

def all_features():
    return [ImageMagick()]
