"""
Histograms
"""

#*****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.plot.primitive import GraphicPrimitive
from sage.plot.plot import minmax_data, Graphics
from sage.plot.misc import options, rename_keyword

class Histogram(GraphicPrimitive):
    """
    Graphics primitive that represents a histogram.  This takes
    quite a few options as well.

    EXAMPLES::

        sage: from sage.plot.histogram import Histogram
        sage: g = Histogram([1,3,2,0], {}); g
        Histogram defined by a data list of size 4
        sage: type(g)
        <class 'sage.plot.histogram.Histogram'>
        sage: opts = { 'bins':20, 'label':'mydata'}
        sage: g = Histogram([random() for _ in range(500)], opts); g
        Histogram defined by a data list of size 500

    We can accept multiple sets of the same length::

        sage: g = Histogram([[1,3,2,0], [4,4,3,3]], {}); g
        Histogram defined by 2 data lists
    """
    def __init__(self, datalist, options):
        """
        Initialize a ``Histogram`` primitive along with
        its options.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: Histogram([10,3,5], {'width':0.7})
            Histogram defined by a data list of size 3
        """
        import numpy as np
        self.datalist=np.asarray(datalist,dtype=float)
        if 'linestyle' in options:
            from sage.plot.misc import get_matplotlib_linestyle
            options['linestyle'] = get_matplotlib_linestyle(
                    options['linestyle'], return_type='long')
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Get minimum and maximum horizontal and vertical ranges
        for the Histogram object.

        EXAMPLES::

            sage: H = histogram([10,3,5], normed=True); h = H[0]
            sage: h.get_minmax_data()
            {'xmax': 10.0, 'xmin': 3.0, 'ymax': 0.4761904761904765, 'ymin': 0}
            sage: G = histogram([random() for _ in range(500)]); g = G[0]
            sage: g.get_minmax_data() # random output
            {'xmax': 0.99729312925213209, 'xmin': 0.00013024562219410285, 'ymax': 61, 'ymin': 0}
            sage: Y = histogram([random()*10 for _ in range(500)], range=[2,8]); y = Y[0]
            sage: ymm = y.get_minmax_data(); ymm['xmax'], ymm['xmin']
            (8.0, 2.0)
            sage: Z = histogram([[1,3,2,0], [4,4,3,3]]); z = Z[0]
            sage: z.get_minmax_data()
            {'xmax': 4.0, 'xmin': 0, 'ymax': 2, 'ymin': 0}
        """
        import numpy
        options=self.options()
        opt=dict(range = options.pop('range',None),
                 bins = options.pop('bins',None),
                 normed = options.pop('normed',None),
                 weights = options.pop('weights', None))
 
        #check to see if a list of datasets
        if not hasattr(self.datalist[0],'__contains__' ):
            ydata,xdata=numpy.histogram(self.datalist, **opt)
            return minmax_data(xdata,[0]+list(ydata), dict=True)
        else:
            m = { 'xmax': 0, 'xmin':0, 'ymax':0, 'ymin':0}
            if not options.pop('stacked',None):
                for d in self.datalist:
                    ydata, xdata = numpy.histogram(d,**opt)
                    m['xmax'] = max([m['xmax']] + list(xdata))
                    m['xmin'] = min([m['xmin']] + list(xdata))
                    m['ymax'] = max([m['ymax']] + list(ydata))
                return m
            else:
                for d in self.datalist:
                    ydata, xdata = numpy.histogram(d,**opt)
                    m['xmax'] = max([m['xmax']] + list(xdata))
                    m['xmin'] = min([m['xmin']] + list(xdata))
                    m['ymax'] = m['ymax'] + max(list(ydata))
                return m
     
    def _allowed_options(self):
        """
        Return the allowed options with descriptions for this graphics
        primitive. This is used in displaying an error message when the
        user gives an option that doesn't make sense.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: g = Histogram( [1,3,2,0], {})
            sage: L = list(sorted(g._allowed_options().items()))
            sage: L[0]
            ('align',
             'How the bars align inside of each bin. Acceptable values are "left", "right" or "mid".')
            sage: L[-1]
            ('zorder', 'The layer level to draw the histogram')
            sage: len(L)
            18
        """
        return {'color': 'The color of the face of the bars or list of colors if multiple data sets are given.',
                'edgecolor':'The color of the border of each bar.',
                'alpha': 'How transparent the plot is',
                'hue':'The color of the bars given as a hue.',
                'fill':'(True or False, default True) Whether to fill the bars',
                'hatch': 'What symbol to fill with - one of "/", "\\", "|", "-", "+", "x", "o", "O", ".", "*"',
                'linewidth':'Width of the lines defining the bars',
                'linestyle': "One of 'solid' or '-', 'dashed' or '--', 'dotted' or ':', 'dashdot' or '-.'",
                'zorder':'The layer level to draw the histogram',
                'bins': 'The number of sections in which to divide the range. Also can be a sequence of points within the range that create the partition.', 
                'align': 'How the bars align inside of each bin. Acceptable values are "left", "right" or "mid".',
                'rwidth': 'The relative width of the bars as a fraction of the bin width',
                'cumulative': '(True or False) If True, then a histogram is computed in which each bin gives the counts in that bin plus all bins for smaller values.  Negative values give a reversed direction of accumulation.',
                'range': 'A list [min, max] which define the range of the histogram. Values outside of this range are treated as outliers and omitted from counts.', 
                'normed': '(True or False) If True, the counts are normalized to form a probability density. (n/(len(x)*dbin)', 
                'weights': 'A sequence of weights the same length as the data list. If supplied, then each value contributes its associated weight to the bin count.',
                'stacked': '(True or False) If True, multiple data are stacked on top of each other.',
                'label': 'A string label for each data list given.'}

    def _repr_(self):
        """
        Return text representation of this histogram graphics primitive.

        EXAMPLES::

            sage: from sage.plot.histogram import Histogram
            sage: g = Histogram( [1,3,2,0], {})
            sage: g._repr_()
            'Histogram defined by a data list of size 4'
            sage: g = Histogram( [[1,1,2,3], [1,3,2,0]], {})
            sage: g._repr_()
            'Histogram defined by 2 data lists'
        """
        L = len(self.datalist)
        if not hasattr(self.datalist[0],'__contains__' ):
            return "Histogram defined by a data list of size {}".format(L)
        else:
            return "Histogram defined by {} data lists".format(L)

    def _render_on_subplot(self, subplot):
        """
        Render this bar chart graphics primitive on a matplotlib subplot
        object.

        EXAMPLES:

        This rendering happens implicitly when the following command
        is executed::

            sage: histogram([1,2,10])  # indirect doctest
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        #check to see if a list of datasets
        if not hasattr(self.datalist[0],'__contains__' ):
            subplot.hist(self.datalist, **options)
        else:
            subplot.hist(self.datalist.transpose(), **options)


@options(aspect_ratio='automatic',align='mid', weights=None, range=None, bins=10, normed=False, edgecolor='black')
def histogram(datalist, **options):
    """
    Computes and draws the histogram for list(s) of numerical data.
    See examples for the many options; even more customization is
    available using matplotlib directly.

    INPUT:

    - ``datalist`` -- A list, or a list of lists, of numerical data
    - ``align`` -- (default: "mid") How the bars align inside of each bin.
      Acceptable values are "left", "right" or "mid"
    - ``alpha`` -- (float in [0,1], default: 1) The transparency of the plot
    - ``bins`` -- The number of sections in which to divide the range. Also
      can be a sequence of points within the range that create the
      partition
    - ``color`` -- The color of the face of the bars or list of colors if
      multiple data sets are given
    - ``cumulative`` -- (boolean - default: False) If True, then
      a histogram is computed in which each bin gives the counts in that
      bin plus all bins for smaller values.  Negative values give
      a reversed direction of accumulation
    - ``edgecolor`` -- The color of the border of each bar
    - ``fill`` -- (boolean - default: True) Whether to fill the bars
    - ``hatch`` -- (default: None) symbol to fill the bars with - one of
      "/", "\\", "|", "-", "+", "x", "o", "O", ".", "*", "" (or None)
    - ``hue`` -- The color of the bars given as a hue. See
      :mod:`~sage.plot.colors.hue` for more information on the hue
    - ``label`` -- A string label for each data list given
    - ``linewidth`` -- (float) width of the lines defining the bars
    - ``linestyle`` -- (default: 'solid') Style of the line. One of 'solid'
      or '-', 'dashed' or '--', 'dotted' or ':', 'dashdot' or '-.'
    - ``normed`` -- (boolean - default: False) If True, the counts are
      normalized to form a probability density.
    - ``range`` -- A list [min, max] which define the range of the
      histogram. Values outside of this range are treated as outliers and
      omitted from counts
    - ``rwidth`` -- (float in [0,1], default: 1) The relative width of the bars
      as a fraction of the bin width
    - ``stacked`` -- (boolean - default: False) If True, multiple data are
      stacked on top of each other
    - ``weights`` -- (list) A sequence of weights the same length as the data
      list. If supplied, then each value contributes its associated weight
      to the bin count
    - ``zorder`` -- (integer) the layer level at which to draw the histogram

    .. NOTE::

        The ``weights`` option works only with a single list. List of lists
        representing multiple data are not supported.

    EXAMPLES:

    A very basic histogram for four data points::

        sage: histogram([1,2,3,4], bins=2)
        Graphics object consisting of 1 graphics primitive

    We can see how the histogram compares to various distributions.
    Note the use of the ``normed`` keyword to guarantee the plot
    looks like the probability density function::

        sage: nv = normalvariate
        sage: H = histogram([nv(0,1) for _ in range(1000)], bins=20, normed=True, range=[-5,5])
        sage: P = plot( 1/sqrt(2*pi)*e^(-x^2/2), (x,-5,5), color='red', linestyle='--')
        sage: H+P
        Graphics object consisting of 2 graphics primitives

    There are many options one can use with histograms.  Some of these
    control the presentation of the data, even if it is boring::

        sage: histogram(list(range(100)), color=(1,0,0), label='mydata',\
              rwidth=.5, align="right")
        Graphics object consisting of 1 graphics primitive

    This includes many usual matplotlib styling options::

        sage: T = RealDistribution('lognormal', [0,1])
        sage: histogram( [T.get_random_element() for _ in range(100)], alpha=0.3,\
              edgecolor='red', fill=False, linestyle='dashed', hatch='O', linewidth=5)
        Graphics object consisting of 1 graphics primitive
        sage: histogram( [T.get_random_element() for _ in range(100)],linestyle='-.')
        Graphics object consisting of 1 graphics primitive

    We can do several data sets at once if desired::

        sage: histogram([srange(0,1,.1)*10, [nv(0, 1) for _ in range(100)]], color=['red','green'], bins=5)
        Graphics object consisting of 1 graphics primitive

    We have the option of stacking the data sets too::

        sage: histogram([ [1,1,1,1,2,2,2,3,3,3], [4,4,4,4,3,3,3,2,2,2] ], stacked=True, color=['blue', 'red'])
        Graphics object consisting of 1 graphics primitive

    It is possible to use weights with the histogram as well::

        sage: histogram(list(range(10)), bins=3, weights=[1,2,3,4,5,5,4,3,2,1])
        Graphics object consisting of 1 graphics primitive
    """
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Histogram(datalist, options=options))
    return g
