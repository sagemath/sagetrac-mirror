"""
Wraps some Python Imaging Library functions in Sage.

Uses PIL version >= 1.1.6 and matplotlib. (This also assumes some system
libraries are installed, such as libjpeg, and will not function properly
otherwise.)

This module contains various PIL functions adapted from the
PIL handbook (http://www.pythonware.com/library/pil/handbook/index.htm,
whose license is stated below) and other sources.

############################################################################
The Python Imaging Library is:

   Copyright (c) 1997-2005 by Secret Labs AB
   Copyright (c) 1995-2005 by Fredrik Lundh

By obtaining, using, and/or copying this software and/or its associated
documentation, you agree that you have read, understood, and will comply
with the following terms and conditions:

Permission to use, copy, modify, and distribute this software and its
associated documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appears in all copies, and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Secret Labs AB or the author not be used
in advertising or publicity pertaining to distribution of the software without
specific, written prior permission.

SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR BE LIABLE FOR ANY SPECIAL,
INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

http://www.pythonware.com/library/pil/handbook/license.htm
############################################################################


AUTHORS:

- 2009-8:  David Joyner <wdjoyner@gmail.com> - original module

.. TODO::

    - Provide an interface with Graphics objects
      (perhaps one could enter such an object instead of a
      string representing a path to a graphics file in a
      suitable format)?
    - Set up a "wrapper class" for PIL (say, "class PILimage(SageObject):"?)
      and move the wrapper functions below to methods of it?

The test images are taken by David Joyner and are placed in
the public domain (pil-test-image1.png, pil-test-image2.png).
"""
################################################################################
#       Copyright (C) 2005 Secret Labs AB and Fredrik Lundh
#                     2009 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL), version 2
#  or later (at your option). The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
################################################################################

import os
import Image
import ImageFont
import ImageDraw
import ImageFilter
import urllib


def image_convert(files, format="jpg"):
    """
    Converts to various image formats using PIL.

    Saves result in SAGE_TMP.

    .. NOTE::

        Converting to/from jpg requires, in addition, libjpeg.

    EXAMPLES::

        sage: from sage.media.pil import image_convert

        sage: infile1 = "http://www.sagemath.org/old/art/.html/sage-car.png"
        sage: image_convert(infile1)  # optional -- internet

        sage: infile2 = "http://www.sagemath.org/old/art/.html/sage-icon7.png"

        sage: files = [infile1, infile2]
        sage: image_convert(files,"gif")  # optional -- internet

    From http://www.pythonware.com/library/pil/handbook/introduction.htm
    """
    if not(type(files) == list):  # assume files is a single file
        files = [files]

    for infile in files:
        path, filename = os.path.split(infile)
        pre, ext = os.path.splitext(filename)
        outfile = os.environ['SAGE_ROOT'] + "/tmp/" + pre + "." + format
        if infile != outfile:
            try:
                Image.open(infile).save(outfile)
            except IOError:
                #print "cannot convert", infile
                pass


def image_blend(infile1, infile2, outfile, mask_ratio=0.5):
    """
    Wraps PIL's blend function (which I cannot find in the manual).

    infile1, infile2 - should be of the same mode and
       the same size. However, cropping (starting at the upper left-hand
       corner) and converting is done just in case they are different.
    outfile - a string representing a legal file name of
       the correct format in a directory you have write access to.
    mask_ration - the transparency layer for the mask image infile1.

    EXAMPLES::

        sage: from sage.media.pil import image_blend

        sage: infile1 = "http://www.sagemath.org/old/art/.html/sage-car.png"
        sage: infile2 = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.png"
        sage: image_blend(infile1, infile2, outfile)
    """
    mask = file_loader(infile1)
    img = file_loader(infile2)
    xsize1, ysize1 = mask.size
    xsize2, ysize2 = img.size
    xsize = min(xsize1, xsize2)
    ysize = min(ysize1, ysize2)
    im1 = mask.crop((0, 0, xsize, ysize))
    im2 = img.crop((0, 0, xsize, ysize))
    im1 = im1.convert("RGBA")
    im2 = im2.convert("RGBA")
    img = Image.blend(im1, im2, 0.5)
    img.save(outfile)


def image_roll(infile, outfile, delta):
    """
    'Roll' an image sideways

    EXAMPLES::

        sage: from sage.media.pil import image_roll

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.png"
        sage: image_roll(infile, outfile, 300)
    """
    im = file_loader(infile)
    xsize, ysize = im.size
    delta = delta % xsize
    if delta == 0:
        return im
    part1 = im.crop((0, 0, delta, ysize))
    part2 = im.crop((delta, 0, xsize, ysize))
    im.paste(part2, (0, 0, xsize - delta, ysize))
    im.paste(part1, (xsize - delta, 0, xsize, ysize))
    im.save(outfile)


def image_crop(infile, outfile, xmin, ymin, xmax, ymax):
    """
    Crops the image along the rectangle with lower left-hand corner
    at (xmin, ymin) and upper right-hand corner at (xmax, ymax).

    Wraps the PIL crop function.

    EXAMPLES::

        sage: from sage.media.pil import image_crop

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: image_id(infile)[0][1:]
        ['PNG', '891x837', 'RGB']
        sage: xmin, ymin, xmax, ymax = 10,10, 1000, 700
        sage: outfile = SAGE_TMP+"/pil-test-image1a.png"
        sage: image_crop(infile, outfile, xmin, ymin, xmax, ymax)
    """
    im = file_loader(infile)
    part = im.crop((xmin, ymin, xmax, ymax))
    part.save(outfile)


def image_id(files):
    """
    Returns the file name, graphics format, size, and color "mode".

    From http://www.pythonware.com/library/pil/handbook/introduction.htm

    EXAMPLES::

        sage: from sage.media.pil import image_id

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: image_id(infile)[0][1:]
        ['PNG', '891x837', 'RGB']

        sage: infile2 =  SAGE_ROOT+'/local/lib/python2.7/site-packages/matplotlib/mpl-data/sample_data/lena.jpg'
        sage: image_id(infile2)[0][1:]
        ['JPEG', '512x512', 'RGB']
    """
    if not(type(files) == list):  # assume files is a single file
        files = [files]
    out = []
    for infile in files:
        im = file_loader(infile)
        out += [[infile, im.format, "%dx%d" % im.size, im.mode]]
    return out


def image_filter(infile, outfile, filter_type="BLUR"):
    """
    The current version of the library provides the following set of predefined
    image enhancement filters:

    BLUR, CONTOUR, DETAIL, EDGE_ENHANCE, EDGE_ENHANCE_MORE, EMBOSS, FIND_EDGES,
    SMOOTH, SMOOTH_MORE, and SHARPEN.

    The string X must take on of these values above.

    EXAMPLES::

        sage: from sage.media.pil import image_filter

        sage: infile1 = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: image_convert(infile1, format = "jpg")
        sage: infile2 = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.jpg"
        sage: image_filter(infile2, outfile, filter_type="EDGE_ENHANCE_MORE")
        sage: outfile = SAGE_TMP+"/pil-test-image1b.jpg"
        sage: image_filter(infile2, outfile, filter_type="BLUR")
        sage: outfile = SAGE_TMP+"/pil-test-image1c.jpg"
        sage: image_filter(infile2, outfile, filter_type="SHARPEN")
    """
    im = file_loader(infile)
    if filter_type == "BLUR":
        im1 = im.filter(ImageFilter.BLUR)
    if filter_type == "CONTOUR":
        im1 = im.filter(ImageFilter.CONTOUR)
    if filter_type == "DETAIL":
        im1 = im.filter(ImageFilter.DETAIL)
    if filter_type == "EDGE_ENHANCE":
        im1 = im.filter(ImageFilter.EDGE_ENHANCE)
    if filter_type == "EDGE_ENHANCE_MORE":
        im1 = im.filter(ImageFilter.EDGE_ENHANCE_MORE)
    if filter_type == "EMBOSS":
        im1 = im.filter(ImageFilter.EMBOSS)
    if filter_type == "FIND_EDGES":
        im1 = im.filter(ImageFilter.FIND_EDGES)
    if filter_type == "SMOOTH":
        im1 = im.filter(ImageFilter.SMOOTH)
    if filter_type == "SMOOTH_MORE":
        im1 = im.filter(ImageFilter.SMOOTH_MORE)
    if filter_type == "SHARPEN":
        im1 = im.filter(ImageFilter.SHARPEN)
    im1.save(outfile)


def image_filter_median(infile, outfile, size=3):
    """
    For each pixel in the input image, this filter copies the median
    pixel value from a (size, size) environment to the output image.
    It appears that size must be an odd integer > 0.

    Other filters like this one are MinFilter(size=3), MaxFilter(size=3),
    and ModeFilter(size=3).

    EXAMPLES::

        sage: from sage.media.pil import image_filter_median

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.png"
        sage: image_filter_median(infile, outfile)

    REFERENCES:

    .. http://www.pythonware.com/library/pil/handbook/imagefilter.htm
    """
    im = file_loader(infile)
    im1 = im.filter(ImageFilter.MedianFilter(size))
    im1.save(outfile)


def image_draw_text(infile, outfile, pt=(10, 10),
                    text=None, font_type="arial.pil", font_size=18):
    """
    A bitmap font is represented by a "pil" extension. A truetype font is
    represented by a "ttf" extension, such as
    font=ImageFont.truetype("symbol.ttf", 16, encoding="symb").

    EXAMPLES::

        sage: from sage.media.pil import image_draw_text
        sage: from matplotlib import font_manager

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.png"
        sage: fontstr = font_manager.findfont("courier") # truetype font "close" to Courier
        sage: vera = ImageFont.truetype(fontstr, 18)
        sage: image_draw_text(infile, outfile, (10,10), "Annapolis Mall ceiling", fontstr)

    You should see "Annapolis Mall ceiling" in the outfile image. (Moreover,
    the value of the variable fontstr is equal to the string
    SAGE_ROOT+"/local/lib/python/site-packages/matplotlib/mpl-data/fonts/ttf/Vera.ttf".)

    REFERENCES:

    .. http://www.pythonware.com/library/pil/handbook/imagefilter.htm
    """
    im = file_loader(infile)
    draw = ImageDraw.Draw(im)
    #print font_type,  font_type[int(-3):], font_type[int(-3):] == "ttf"
    if font_type[int(-3):] == "ttf":
        fonttype = ImageFont.truetype(font_type, font_size)
    if font_type[int(-3):] == "pil":
        fonttype = ImageFont.load(font_type)
    draw.text(pt, text, font=fonttype)
    im.save(outfile)


def image_pixel_color(infile, x, y):
    """
    Returns the color tuple of the image infile at the (x,y)-th pixel,
    where (0,0) is the upper left-hand pixel.

    EXAMPLES::

        sage: from sage.media.pil import image_pixel_color

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: image_pixel_color(infile, 1, 1)
        (255, 255, 255)

    Here the 4th coordinate represents the alpha channel since the image
    has RGBA mode.

    TESTS::

        sage: image_pixel_color(infile, 6913, 6914)
        Traceback (most recent call last):
        ...
        ValueError: The pixel values are off the 891 x 837 image
    """
    im = file_loader(infile)
    width, height = im.size
    if not(0 <= x and x <= width and 0 <= y and y <= height):
        raise ValueError("The pixel values are off the %s x %s image" % (width, height))
    return im.getpixel((int(x), int(y)))


def image_statistics(infile):
    """
    From http://www.pythonware.com/library/pil/handbook/imagestat.htm
    Returns the pixel count, mean pixel level, median pixel level,
    and its standard deviation.

    EXAMPLES::

        sage: from sage.media.pil import image_statistics

        sage: infile = "http://www.sagemath.org/old/art/.html/sage-icon7.png"
        sage: image_statistics(infile)
        [[745767, 745767, 745767], [202.52117350325236, 213.39012184770846, 206.43753477962957], [254, 254, 254], [78.80113084931638, 62.11183125319092, 75.20923237875257]]

        sage: infile2 = SAGE_ROOT+'/devel/sage-main/doc/fr/a_tour_of_sage/eigen_plot.png'
        sage: image_statistics(infile2)
        [[300000, 300000, 300000, 300000], [250.5215, 250.5215, 253.53336333333334, 255.0], [255, 255, 255, 255], [31.64086257384064, 31.64086257384064, 17.922602495768185, 0.0]]
    """
    from ImageStat import Stat
    im = file_loader(infile)
    stats = Stat(im)
    return [stats.count, stats.mean, stats.median, stats.stddev]


def image_autocontrast(infile, outfile):
    """
    Wraps ImageOps' autocontrast method,
    http://www.pythonware.com/library/pil/handbook/imageops.htm,
    which states:
    "Maximize (normalize) image contrast. This function calculates a
    histogram of the input image, removes cutoff percent of the lightest and
    darkest pixels from the histogram, and remaps the image so that the darkest
    remaining pixel becomes black (0), and the lightest becomes white (255)."

    EXAMPLES::

        sage: from sage.media.pil import image_autocontrast

        sage: infile1 = "http://sage.math.washington.edu/home/wdj/art/thumbs/heart.jpg"
        sage: outfile = SAGE_TMP+"/pil-test-image1a.jpg"
        sage: image_autocontrast(infile1, outfile)    # optional -- internet

        sage: infile2 =  SAGE_ROOT+'/local/lib/python2.7/site-packages/matplotlib/mpl-data/sample_data/lena.jpg'
        sage: outfile = SAGE_TMP+"/pil-test-image1a.jpg"
        sage: image_autocontrast(infile2, outfile)
    """
    from ImageOps import autocontrast
    im = file_loader(infile)
    autocontrast(im).save(outfile)


def file_loader(infile):
    """
    Import a file, either local or from the web

    EXAMPLES::

        sage: from sage.media.pil import file_loader

        sage: infile1 = "http://www.sagemath.org/old/art/.html/sage-car.png"
        sage: car = file_loader(infile1)  # optional -- internet
        sage: add(car.histogram())  # optional -- internet
        1505496

        sage: infile2 = SAGE_ROOT+'/devel/sage-main/doc/fr/a_tour_of_sage/eigen_plot.png'
        sage: im = file_loader(infile2)
        sage: add(im.histogram())
        1200000
    """
    try:
        if infile[:4] == "http":
            im1 = urllib.URLopener()
            im2 = im1.retrieve(infile)[0]
            im = Image.open(im2)
        else:
            im = Image.open(infile)
        return im
    except IOError:
        raise RuntimeError('cannot read the file')
