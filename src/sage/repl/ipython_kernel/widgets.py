r"""
Widgets to be used for the Sage Jupyter notebook

These are all based on widgets from ``ipywidgets``, changing or
combining existing widgets.

TESTS:

We need to setup a proper test environment for widgets::

    sage: from ipywidgets.widgets.tests import setup_test_comm
    sage: setup_test_comm()
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import unicode_literals

from ipywidgets.widgets import (IntSlider, IntRangeSlider, FloatSlider,
        FloatRangeSlider, Text, Textarea, ColorPicker, HTML,
        Label, HBox, VBox, ValueWidget)
from traitlets import List, Unicode, link

from sage.misc.sage_eval import sage_eval
from sage.repl.user_globals import get_globals
from sage.plot.colors import Color


class HTMLText(HTML):
    """
    A HTML widget whose ``description`` is always empty.

    This is used to display arbitrary HTML text in interacts without
    a label. The :func:`text_control` function from SageNB is an alias
    of :class:`HTMLText`.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import HTMLText
        sage: w = HTMLText("Hello")
        sage: w.description
        u''
        sage: w.description = "text"
        sage: w.description
        u''
    """
    @property
    def description(self):
        """
        Always return empty string.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.widgets import HTMLText
            sage: w = HTMLText("Hello")
            sage: w.description
            u''
        """
        return ''

    @description.setter
    def description(self, value):
        """
        Do not set anything.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.widgets import HTMLText
            sage: w = HTMLText("Hello")
            sage: w.description = "text"
            sage: w.description
            u''
        """
        pass

class TransformWidget(object):
    """
    A mixin class for a widget to transform the bare widget value for
    use in interactive functions.

    INPUT:

    - ``transform`` -- a one-argument function which transforms the
      value of the widget for use by an interactive function.

    - other arguments are passed to the base class

    EXAMPLES::

        sage: from ipywidgets import ToggleButtons
        sage: from sage.repl.ipython_kernel.widgets import TransformWidget
        sage: class TransformToggleButtons(TransformWidget, ToggleButtons): pass
        sage: w = TransformToggleButtons(options=["pi", "e"], transform=lambda x: x+x)
        sage: w
        TransformToggleButtons(value='pi', options=['pi', 'e'])
        sage: w.get_interact_value()
        'pipi'
    """
    def __init__(self, *args, **kwds):
        """
        Construct a :class:`TransformWidget`.

        TESTS::

            sage: from sage.repl.ipython_kernel.widgets import TransformWidget
            sage: w = TransformWidget(transform=dict)
            sage: w._TransformWidget__transform
            <type 'dict'>
        """
        self.__transform = kwds.pop("transform", None)
        return super(TransformWidget, self).__init__(*args, **kwds)

    def get_value(self):
        """
        Return ``self.value``.

        This is meant to be overridden by sub-classes to change the
        input of the transform function.

        EXAMPLES::

            sage: from ipywidgets import ColorPicker
            sage: from sage.repl.ipython_kernel.widgets import TransformWidget
            sage: class TransformColorPicker(TransformWidget, ColorPicker): pass
            sage: TransformColorPicker(value="red").get_value()
            'red'
        """
        return self.value

    def get_interact_value(self):
        """
        Return the transformed value of this widget, by calling
        the ``transform`` function.

        EXAMPLES:

            sage: from ipywidgets import Checkbox
            sage: from sage.repl.ipython_kernel.widgets import TransformWidget
            sage: class TransformCheckbox(TransformWidget, Checkbox): pass
            sage: w = TransformCheckbox(value=True, transform=int); w
            TransformCheckbox(value=True)
            sage: w.get_interact_value()
            1
        """
        v = self.get_value()
        f = self.__transform
        if f is None:
            return v
        else:
            return f(v)


class EvalWidget(TransformWidget):
    """
    A mixin class for a widget to evaluate (using :func:`sage_eval`) the
    widget value and possibly transform it like :class:`TransformWidget`.

    EXAMPLES::

        sage: from ipywidgets import ToggleButtons
        sage: from sage.repl.ipython_kernel.widgets import EvalWidget
        sage: class EvalToggleButtons(EvalWidget, ToggleButtons): pass
        sage: w = EvalToggleButtons(options=["pi", "e"], transform=lambda x: x+x)
        sage: w
        EvalToggleButtons(value='pi', options=['pi', 'e'])
        sage: w.get_interact_value()
        2*pi
    """
    def get_value(self):
        """
        Evaluate the bare widget value using :func:`sage_eval`.

        EXAMPLES:

            sage: from ipywidgets import Dropdown
            sage: from sage.repl.ipython_kernel.widgets import EvalWidget
            sage: class EvalDropdown(EvalWidget, Dropdown): pass
            sage: w = EvalDropdown(options=["the_answer"], transform=RR)
            sage: w
            EvalDropdown(value='the_answer', options=['the_answer'])
            sage: the_answer = 42
            sage: w.get_value()
            42
            sage: w.get_interact_value()
            42.0000000000000
        """
        return sage_eval(self.value, get_globals())


class TransformIntSlider(TransformWidget, IntSlider):
    """
    An :class:`ipywidgets.IntSlider` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformIntSlider
        sage: w = TransformIntSlider(min=0, max=100, value=7, transform=lambda x: x^2)
        sage: w
        TransformIntSlider(value=7, min=0, max=100, step=1)
        sage: w.get_interact_value()
        49
    """
    pass


class TransformFloatSlider(TransformWidget, FloatSlider):
    """
    A :class:`ipywidgets.FloatSlider` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformFloatSlider
        sage: w = TransformFloatSlider(min=0, max=100, value=7, transform=lambda x: sqrt(x))
        sage: w
        TransformFloatSlider(value=7.0, min=0.0, max=100.0, step=0.1)
        sage: w.get_interact_value()
        2.6457513110645907
    """
    pass


class TransformIntRangeSlider(TransformWidget, IntRangeSlider):
    """
    An :class:`ipywidgets.IntRangeSlider` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformIntRangeSlider
        sage: w = TransformIntRangeSlider(min=0, max=100, value=(7,9), transform=lambda x: x[1]-x[0])
        sage: w
        TransformIntRangeSlider(value=(7, 9), min=0, max=100, step=1)
        sage: w.get_interact_value()
        2
    """
    pass


class TransformFloatRangeSlider(TransformWidget, FloatRangeSlider):
    """
    An :class:`ipywidgets.FloatRangeSlider` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformFloatRangeSlider
        sage: w = TransformFloatRangeSlider(min=0, max=100, value=(7,9), transform=lambda x: x[1]-x[0])
        sage: w
        TransformFloatRangeSlider(value=(7.0, 9.0), min=0.0, max=100.0, step=1.0)
        sage: w.get_interact_value()
        2.0
    """
    pass


class TransformText(TransformWidget, Text):
    """
    A :class:`ipywidgets.Text` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformText
        sage: w = TransformText(value="hello", transform=lambda x: x+x)
        sage: w
        TransformText(value=u'hello')
        sage: w.get_interact_value()
        u'hellohello'
    """
    pass


class TransformTextarea(TransformWidget, Textarea):
    """
    A :class:`ipywidgets.Textarea` widget with an optional
    transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import TransformTextarea
        sage: w = TransformTextarea(value="hello", transform=lambda x: x+x)
        sage: w
        TransformTextarea(value=u'hello')
        sage: w.get_interact_value()
        u'hellohello'
    """
    pass


class EvalText(EvalWidget, Text):
    """
    A :class:`ipywidgets.Text` widget which evaluates (using
    :func:`sage_eval`) its contents and applies an optional transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import EvalText
        sage: w = EvalText(value="pi", transform=lambda x: x^2)
        sage: w
        EvalText(value=u'pi')
        sage: w.get_interact_value()
        pi^2
    """
    pass


class EvalTextarea(EvalWidget, Textarea):
    """
    A :class:`ipywidgets.Textarea` widget which evaluates (using
    :func:`sage_eval`) its contents and applies an optional transformation.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import EvalTextarea
        sage: w = EvalTextarea(value="pi", transform=lambda x: x^2)
        sage: w
        EvalTextarea(value=u'pi')
        sage: w.get_interact_value()
        pi^2
    """
    pass


class SageColorPicker(ColorPicker):
    """
    A color picker widget return a Sage :class:`Color`.

    EXAMPLES::

        sage: from sage.repl.ipython_kernel.widgets import SageColorPicker
        sage: SageColorPicker()
        SageColorPicker(value='black')
    """
    def get_interact_value(self):
        """
        Return a Sage :class:`Color` corresponding to the value of this
        widget.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.widgets import SageColorPicker
            sage: SageColorPicker().get_interact_value()
            RGB color (0.0, 0.0, 0.0)
        """
        return Color(self.value)


class Grid(TransformWidget, HBox, ValueWidget):
    """
    A square grid of widgets whose value is a list of lists of the
    values of the individual widgets.

    This is usually created using the :func:`input_grid` function.

    EXAMPLES::

        sage: from ipywidgets import Text
        sage: from sage.repl.ipython_kernel.widgets import Grid
        sage: w = Grid(2, 2, lambda i,j: Text(value="%s,%s"%(i,j)))
        sage: w
        Grid(value=[[u'0,0', u'0,1'], [u'1,0', u'1,1']], children=(Label(value=u''), VBox(children=(Text(value=u'0,0'), Text(value=u'1,0'))), VBox(children=(Text(value=u'0,1'), Text(value=u'1,1')))))
        sage: w.get_interact_value()
        [[u'0,0', u'0,1'], [u'1,0', u'1,1']]
    """
    value = List()
    description = Unicode()

    def __init__(self, nrows, ncols, make_widget, description=u"", transform=None):
        """
        Create a :class:`Grid` widget.

        INPUT:

        - ``nrows``, ``ncols`` -- number of rows and columns in the grid

        - ``make_widget`` -- a function of two arguments ``(i,j)``
          returning the widget to be placed at position ``(i,j)``.

        - ``description`` -- an optional label.

        - ``transform`` -- an optional transformation, see :class:`TransformWidget`.

        EXAMPLES::

            sage: from sage.repl.ipython_kernel.widgets import Grid, EvalText
            sage: w = Grid(2, 2, lambda i,j: EvalText(str(j+4*i)),
            ....:         description="2x2 matrix", transform=matrix)
            sage: w
            Grid(value=[[0, 1], [4, 5]], description=u'2x2 matrix', children=(Label(value=u'2x2 matrix'), VBox(children=(EvalText(value=u'0'), EvalText(value=u'4'))), VBox(children=(EvalText(value=u'1'), EvalText(value=u'5')))))
            sage: w.get_interact_value()
            [0 1]
            [4 5]

        TESTS::

            sage: w = Grid(0, 1, lambda i,j: EvalText())
            Traceback (most recent call last):
            ...
            ValueError: Grid requires a positive number of rows and columns
        """
        if nrows < 1 or ncols < 1:
            raise ValueError("Grid requires a positive number of rows and columns")
        super(Grid, self).__init__(transform=transform)

        label = Label(description)
        link((label, "value"), (self, "description"))

        self.cols = []
        for j in range(ncols):
            col = VBox()
            widgets = []
            for i in range(nrows):
                w = make_widget(i, j)
                w.observe(self._update, names="value")
                widgets.append(w)
            col.children = widgets
            self.cols.append(col)
        self.children = [label] + self.cols
        self._update()

    def _update(self, *args):
        """
        Compute the ``value`` of the grid.

        Thanks to traitlets magic, this is automatically called
        whenever one of the child widgets changes.

        EXAMPLES::

            sage: from ipywidgets import Text
            sage: from sage.repl.ipython_kernel.widgets import Grid
            sage: w = Grid(2, 2, lambda i,j: Text(value="%s,%s"%(i,j)))
            sage: w._update()
            sage: w.value
            [[u'0,0', u'0,1'], [u'1,0', u'1,1']]
            sage: w.cols[0].children[0].value = "abc"
            sage: w.value
            [[u'abc', u'0,1'], [u'1,0', u'1,1']]
        """
        v = [[]]
        for col in self.cols:
            for i in range(len(col.children)):
                if i >= len(v):
                    v.append([])
                v[i].append(col.children[i].get_interact_value())
        self.value = v
