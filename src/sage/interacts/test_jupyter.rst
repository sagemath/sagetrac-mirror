.. -*- coding: utf-8 -*-

.. linkall

Test the Sage interact library in Jupyter

We need to setup a proper test environment for widgets::

    sage: from ipywidgets.widgets.tests import setup_test_comm
    sage: setup_test_comm()

Make sure that we use the Jupyter interacts::

    sage: from sage.repl.ipython_kernel.all_jupyter import *

Function to test library interacts. This will "display" the widgets
textually and it will actually run the interactive function. Graphics
which are produced by the interact are not displayed, but text is. ::

    sage: def test(lib_interact):
    ....:     f = lib_interact.f      # Underlying non-wrapped function
    ....:     w = interact(f).widget  # "interactive" widget
    ....:     kwargs = {widget._kwarg: widget.get_interact_value()
    ....:             for widget in w.kwargs_widgets}
    ....:     return f(**kwargs)

This is just to test that failures in the interact are actually seen::

    sage: from sage.interacts.library import library_interact
    sage: @library_interact
    ....: def failure():
    ....:     raise Exception("gotcha")
    sage: test(failure)
    Traceback (most recent call last):
    ...
    Exception: gotcha

Test all interacts from the Sage interact library::

    sage: test(interacts.algebra.polar_prime_spiral)  # long time
    Interactive function <function polar_prime_spiral at ...> with 6 widgets
      interval: IntRangeSlider(value=(1, 1000), min=1, max=4000, step=10, description=u'range')
      show_factors: Checkbox(value=True, description=u'show_factors')
      highlight_primes: Checkbox(value=True, description=u'highlight_primes')
      show_curves: Checkbox(value=True, description=u'show_curves')
      n: IntSlider(value=89, min=1, max=200, step=1, description=u'number $n$')
      dpi: IntSlider(value=100, min=10, max=300, step=10, description=u'dpi')
    <h2>Polar Prime Spiral</h2>           <div style="white-space: normal;">          For more information about the factors in the spiral, visit           <a href="http://www.dcs.gla.ac.uk/~jhw/spirals/index.html" target="_blank">          Number Spirals by John Williamson</a>.</div>
    <script type="math/tex">n = 89</script>
    Pink Curve:  <script type="math/tex">n^2 + 8.0</script>
    Green Curve: <script type="math/tex">n^2 + n + -1.0</script>

    sage: test(interacts.calculus.taylor_polynomial)
    Interactive function <function taylor_polynomial at ...> with 3 widgets
      title: HTMLText(value=u'<h2>Taylor polynomial</h2>')
      f: EvalText(value=u'e^(-x)*sin(x)', description=u'$f(x)=$')
      order: SelectionSlider(value=1, options=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], description=u'order')
    <script type="math/tex">f(x)\;=\;e^{\left(-x\right)} \sin\left(x\right)</script>
    <script type="math/tex">\hat{f}(x;0)\;=\;x+\mathcal{O}(x^{2})</script>

    sage: test(interacts.calculus.definite_integral)
    Interactive function <function definite_integral at ...> with 6 widgets
      title: HTMLText(value=u'<h2>Definite integral</h2>')
      f: EvalText(value=u'3*x', description=u'$f(x)=$')
      g: EvalText(value=u'x^2', description=u'$g(x)=$')
      interval: IntRangeSlider(value=(0, 3), min=-10, max=10, step=1, description=u'Interval')
      x_range: IntRangeSlider(value=(0, 3), min=-10, max=10, step=1, description=u'plot range (x)')
      selection: Dropdown(value='f and g', options=['f', 'g', 'f and g', 'f - g'], description=u'Select')
    <script type="math/tex">\int_{0.00}^{3.00}(\color{Blue}{f(x)})\,\mathrm{d}x=\int_{0.00}^{3.00}(3 \, x)\,\mathrm{d}x=13.50</script><br/><script type="math/tex">\int_{0.00}^{3.00}(\color{Green}{g(x)})\,\mathrm{d}x=\int_{0.00}^{3.00}(x^{2})\,\mathrm{d}x=9.00</script>

    sage: test(interacts.calculus.function_derivative)
    Interactive function <function function_derivative at ...> with 4 widgets
      title: HTMLText(value=u'<h2>Derivative grapher</h2>')
      function: EvalText(value=u'x^5-3*x^3+1', description=u'Function:')
      x_range: FloatRangeSlider(value=(-2.0, 2.0), min=-15.0, max=15.0, step=0.1, description=u'Range (x)')
      y_range: FloatRangeSlider(value=(-8.0, 6.0), min=-15.0, max=15.0, step=0.1, description=u'Range (y)')
    <center><script type="math/tex">\color{Blue}{f(x) = x^{5} - 3 \, x^{3} + 1}</script></center>
    <center><script type="math/tex">\color{Green}{f'(x) = 5 \, x^{4} - 9 \, x^{2}}</script></center>
    <center><script type="math/tex">\color{Red}{f''(x) = 20 \, x^{3} - 18 \, x}</script></center>

    sage: test(interacts.calculus.difference_quotient)
    Interactive function <function difference_quotient at ...> with 5 widgets
      title: HTMLText(value=u'<h2>Difference quotient</h2>')
      f: EvalText(value=u'sin(x)', description=u'f(x)')
      interval: FloatRangeSlider(value=(0.0, 10.0), min=0.0, max=10.0, step=0.1, description=u'Range')
      a: IntSlider(value=5, min=0, max=10, step=1, description=u'$a$')
      x0: IntSlider(value=2, min=0, max=10, step=1, description=u'$x_0$ (start point)')
    <h2>Difference Quotient</h2>
    <div style="white-space: normal;">         <a href="http://en.wikipedia.org/wiki/Difference_quotient" target="_blank">         Wikipedia article about difference quotient</a></div>
    <h2>Difference Quotient</h2>
    <br><script type="math/tex">\text{Line's equation:}</script>
    <script type="math/tex">y = 1/3*(x - 5)*(sin(5) - sin(2)) + sin(5)</script><br>
    <script type="math/tex">\text{Slope:}</script>
    <script type="math/tex">k = \frac{f(x_0)-f(a)}{x_0-a} = -0.62274</script><br>

    sage: test(interacts.calculus.quadratic_equation)
    Interactive function <function quadratic_equation at ...> with 3 widgets
      A: IntSlider(value=1, min=-7, max=7, step=1, description=u'A')
      B: IntSlider(value=1, min=-7, max=7, step=1, description=u'B')
      C: IntSlider(value=-2, min=-7, max=7, step=1, description=u'C')
    <h2>The Solutions of the Quadratic Equation</h2>
    <script type="math/tex">x^2 + x - 2 = 0</script>
    <script type="math/tex">Ax^2 + Bx + C = 0</script>
    <script type="math/tex">x = \frac{-B\pm\sqrt{B^2-4AC}}{2A} = \frac{-1\pm\sqrt{1^2-4*1*-2}}{2*1} = \frac{-1\pm\sqrt{\color{Green}{9}}}{2} = \begin{cases}1\\-2\end{cases}</script>

    sage: test(interacts.calculus.trigonometric_properties_triangle)
    Interactive function <function trigonometric_properties_triangle at ...> with 3 widgets
      a0: IntSlider(value=30, min=0, max=360, step=1, description=u'A')
      a1: IntSlider(value=180, min=0, max=360, step=1, description=u'B')
      a2: IntSlider(value=300, min=0, max=360, step=1, description=u'C')
    <h2>Trigonometric Properties of a Triangle</h2>
    <script type="math/tex">\angle A = {60.0}^{\circ},</script> <script type="math/tex">\angle B = {45.0}^{\circ},</script> <script type="math/tex">\angle C = {75.0}^{\circ}</script>
    <script type="math/tex">AB = 1.93185165258,</script>  <script type="math/tex">BC = 1.73205080757,</script>  <script type="math/tex">CA = 1.41421356237</script>
    Area of triangle <script type="math/tex">ABC = 1.18301270189</script>

    sage: test(interacts.calculus.secant_method)
    Interactive function <function secant_method at ...> with 5 widgets
      title: HTMLText(value=u'<h2>Secant method for numerical root finding</h2>')
      f: EvalText(value=u'x^2-2', description=u'f(x)')
      interval: IntRangeSlider(value=(0, 4), min=-5, max=5, step=1, description=u'range')
      d: IntSlider(value=3, min=1, max=16, step=1, description=u'10^-d precision')
      maxn: IntSlider(value=10, min=0, max=15, step=1, description=u'max iterations')
    <script type="math/tex">\text{Precision }h = 10^{-d}=10^{-3}=0.00100</script>
    <script type="math/tex">{c = }1.41440380977</script>
    <script type="math/tex">{f(c) = }0.000538137094544</script>
    <script type="math/tex">6 \text{ iterations}</script>

    sage: test(interacts.calculus.newton_method)
    Interactive function <function newton_method at ...> with 7 widgets
      title: HTMLText(value=u'<h2>Newton method</h2>')
      f: EvalText(value=u'x^2 - 2', description=u'f')
      c: IntSlider(value=6, min=-10, max=10, step=1, description=u'Start ($x$)')
      d: IntSlider(value=3, min=1, max=16, step=1, description=u'$10^{-d}$ precision')
      maxn: IntSlider(value=10, min=0, max=15, step=1, description=u'max iterations')
      interval: IntRangeSlider(value=(0, 6), min=-10, max=10, step=1, description=u'Interval')
      list_steps: Checkbox(value=False, description=u'List steps')
    <script type="math/tex">\text{Precision } 2h = 0.001</script>
    <script type="math/tex">{c = }1.4142141576301823</script>
    <script type="math/tex">{f(c) = }1.6836416461 \times 10^{-06}</script>
    <script type="math/tex">6 \text{ iterations}</script>

    sage: test(interacts.calculus.trapezoid_integration)
    Interactive function <function trapezoid_integration at ...> with 7 widgets
      title: HTMLText(value=u'<h2>Trapezoid integration</h2>')
      f: EvalText(value=u'x^2-5*x + 10', description=u'$f(x)=$')
      n: IntSlider(value=5, min=1, max=100, step=1, description=u'# divisions')
      interval_input: ToggleButtons(value='from slider', options=['from slider', 'from keyboard'], description=u'Integration interval')
      interval_s: IntRangeSlider(value=(0, 8), min=-10, max=10, step=1, description=u'slider: ')
      interval_g: Grid(value=[[0, 8]], description=u'keyboard: ', children=(Label(value=u'keyboard: '), VBox(children=(EvalText(value=u'0'),)), VBox(children=(EvalText(value=u'8'),))))
      output_form: ToggleButtons(value='traditional', options=['traditional', 'table', 'none'], description=u'Computations form')
    Function <script type="math/tex">f(x)=x^{2} - 5 \, x + 10</script>
    Integral value to seven decimal places is: <script type="math/tex">\displaystyle\int_{0.00}^{8.00} {f(x) \, \mathrm{d}x} = 90.666667</script>
    <BLANKLINE>
                <div class="math">
                \begin{align*}
                \int_{0.00}^{8.00} {f(x) \, \mathrm{d}x}
                    & \approx \frac {d}{2} \cdot \left[f(x_0) + 2 f(x_{1}) + 2 f(x_{2}) + 2 f(x_{3}) + 2 f(x_{4}) + f(x_{5})\right] \\
                    & = \frac{1.60}{2} \cdot \left[f(0.00) + 2 f(1.60) + 2 f(3.20) + 2 f(4.80) + 2 f(6.40) + f(8.00)\right] \\
                    & = \frac{1.60}{2} \cdot \left[10.00 + 2\cdot 4.56 + 2\cdot 4.24 + 2\cdot 9.04 + 2\cdot 18.96 + 34.00\right] \\
                    & = 94.08000
                \end{align*}
                </div>
    <BLANKLINE>

    sage: test(interacts.calculus.simpson_integration)
    Interactive function <function simpson_integration at ...> with 7 widgets
      title: HTMLText(value=u'<h2>Simpson integration</h2>')
      f: EvalText(value=u'x*sin(x)+x+1', description=u'$f(x)=$')
      n: IntSlider(value=6, min=2, max=100, step=2, description=u'# divisions')
      interval_input: ToggleButtons(value='from slider', options=['from slider', 'from keyboard'], description=u'Integration interval')
      interval_s: IntRangeSlider(value=(0, 10), min=-10, max=10, step=1, description=u'slider: ')
      interval_g: Grid(value=[[0, 10]], description=u'keyboard: ', children=(Label(value=u'keyboard: '), VBox(children=(EvalText(value=u'0'),)), VBox(children=(EvalText(value=u'10'),))))
      output_form: ToggleButtons(value='traditional', options=['traditional', 'table', 'none'], description=u'Computations form')
    Function <script type="math/tex">f(x)=x \sin\left(x\right) + x + 1</script>
    Integral value to seven decimal places is: <script type="math/tex">\displaystyle\int_{0.00}^{10.00} {f(x) \, \mathrm{d}x} = 67.846694</script>
    <BLANKLINE>
            <div class="math">
            \begin{align*}
            \int_{0.00}^{10.00} {f(x) \, \mathrm{d}x}
                & \approx \frac{d}{3} \cdot \left[ f(x_0) + 4 \cdot f(x_{1}) + 2 \cdot f(x_{2}) + 4 \cdot f(x_{3}) + 2 \cdot f(x_{4}) + 4 \cdot f(x_{5}) + f(x_{6})\right] \\
                & = \frac{1.67}{3} \cdot \left[ f(0.00) +  4 \cdot f(1.67) + 2 \cdot f(3.33) + 4 \cdot f(5.00) + 2 \cdot f(6.67) + 4 \cdot f(8.33) + f(10.00)\right] \\
                & = \frac{1.67}{3} \cdot \left[ 1.00 +  4 \cdot 4.33 + 2 \cdot 3.70 + 4 \cdot 1.21 + 2 \cdot 10.16 + 4 \cdot 16.73  + 5.56\right] \\
                & = 68.506699
            \end{align*}
            </div>
    <BLANKLINE>

    sage: test(interacts.calculus.bisection_method)
    Interactive function <function bisection_method at ...> with 5 widgets
      title: HTMLText(value=u'<h2>Bisection method</h2>')
      f: EvalText(value=u'x^2-2', description=u'f(x)')
      interval: IntRangeSlider(value=(0, 4), min=-5, max=5, step=1, description=u'range')
      d: IntSlider(value=3, min=1, max=8, step=1, description=u'$10^{-d}$ precision')
      maxn: IntSlider(value=10, min=0, max=50, step=1, description=u'max iterations')
    <script type="math/tex">\text{Precision }h = 10^{-d}=10^{-3}=0.00100</script>
    <script type="math/tex">{c = }1.4140625</script>
    <script type="math/tex">{f(c) = }-0.00042724609375</script>
    <script type="math/tex">9 \text{ iterations}</script>

    sage: test(interacts.calculus.riemann_sum)
    Manual interactive function <function riemann_sum at ...> with 9 widgets
      title: HTMLText(value=u'<h2>Riemann integral with random sampling</h2>')
      f: EvalText(value=u'x^2+1', description=u'$f(x)=$')
      n: IntSlider(value=5, min=1, max=30, step=1, description=u'# divisions')
      hr1: HTMLText(value=u'<hr>')
      interval_input: ToggleButtons(value='from slider', options=['from slider', 'from keyboard'], description=u'Integration interval')
      interval_s: IntRangeSlider(value=(0, 2), min=-5, max=10, step=1, description=u'slider: ')
      interval_g: Grid(value=[[0, 2]], description=u'keyboard: ', children=(Label(value=u'keyboard: '), VBox(children=(EvalText(value=u'0'),)), VBox(children=(EvalText(value=u'2'),))))
      hr2: HTMLText(value=u'<hr>')
      list_table: Checkbox(value=False, description=u'List table')
    <small>Adjust your data and click Update button. Click repeatedly for another random values.</small>
    Riemann sum: <script type="math/tex">\displaystyle\sum_{i=1}^{5} f(\eta_i)(x_i-x_{i-1})=...</script>
    Exact value of the integral <script type="math/tex">\displaystyle\int_{0}^{2}x^{2} + 1\,\mathrm{d}x=4.66666666667</script>

    sage: test(interacts.calculus.function_tool)
    Interactive function <function function_tool at ...> with 7 widgets
      f: EvalText(value=u'sin(x)', description=u'f')
      g: EvalText(value=u'cos(x)', description=u'g')
      xrange: IntRangeSlider(value=(0, 1), min=-3, max=3, step=1, description=u'x-range')
      yrange: Text(value=u'auto', description=u'yrange')
      a: IntSlider(value=1, min=-1, max=3, step=1, description=u'a')
      action: ToggleButtons(value='f', options=['f', 'df/dx', 'int f', 'num f', 'den f', '1/f', 'finv', 'f+a', 'f-a', 'f*a', 'f/a', 'f^a', 'f(x+a)', 'f(x*a)', 'f+g', 'f-g', 'f*g', 'f/g', 'f(g)'], description=u'h = ')
      do_plot: Checkbox(value=True, description=u'do_plot')
    <center><font color="red"><script type="math/tex">f = \sin\left(x\right)</script></font></center>
    <center><font color="green"><script type="math/tex">g = \cos\left(x\right)</script></font></center>
    <center><font color="blue"><b><script type="math/tex">h = f = \sin\left(x\right)</script></b></font></center>

    sage: test(interacts.fractals.mandelbrot)
    Interactive function <function mandelbrot at ...> with 6 widgets
      expo: FloatSlider(value=2.0, min=-10.0, max=10.0, step=0.1, description=u'expo')
      iterations: IntSlider(value=20, min=1, max=100, step=1, description=u'# iterations')
      zoom_x: FloatRangeSlider(value=(-2.0, 1.0), min=-2.0, max=2.0, step=0.01, description=u'Zoom X')
      zoom_y: FloatRangeSlider(value=(-1.5, 1.5), min=-2.0, max=2.0, step=0.01, description=u'Zoom Y')
      plot_points: IntSlider(value=150, min=20, max=400, step=20, description=u'plot points')
      dpi: IntSlider(value=80, min=20, max=200, step=10, description=u'dpi')
    <h2>Mandelbrot Fractal</h2>
    Recursive Formula: <script type="math/tex">z \leftarrow z^{2.00} + c</script> for <script type="math/tex">c \in \mathbb{C}</script>

    sage: test(interacts.fractals.julia)
    Interactive function <function julia at ...> with 8 widgets
      expo: FloatSlider(value=2.0, min=-10.0, max=10.0, step=0.1, description=u'expo')
      c_real: FloatSlider(value=0.5, min=-2.0, max=2.0, step=0.01, description=u'real part const.')
      c_imag: FloatSlider(value=0.5, min=-2.0, max=2.0, step=0.01, description=u'imag part const.')
      iterations: IntSlider(value=20, min=1, max=100, step=1, description=u'# iterations')
      zoom_x: FloatRangeSlider(value=(-1.5, 1.5), min=-2.0, max=2.0, step=0.01, description=u'Zoom X')
      zoom_y: FloatRangeSlider(value=(-1.5, 1.5), min=-2.0, max=2.0, step=0.01, description=u'Zoom Y')
      plot_points: IntSlider(value=150, min=20, max=400, step=20, description=u'plot points')
      dpi: IntSlider(value=80, min=20, max=200, step=10, description=u'dpi')
    <h2>Julia Fractal</h2>
    Recursive Formula: <script type="math/tex">z \leftarrow z^{2.00} + (0.50+0.50*\mathbb{I})</script>

    sage: test(interacts.fractals.cellular_automaton)
    Interactive function <function cellular_automaton at ...> with 3 widgets
      N: IntSlider(value=100, min=1, max=500, step=1, description=u'Number of iterations')
      rule_number: IntSlider(value=110, min=0, max=255, step=1, description=u'Rule number')
      size: IntSlider(value=6, min=1, max=11, step=1, description=u'size of graphic')
    <h2>Cellular Automaton</h2><div style="white-space: normal;">"A cellular automaton is a collection of "colored" cells          on a grid of specified shape that evolves through a number of          discrete time steps according to a set of rules based on the          states of neighboring cells." &mdash;          <a target="_blank" href="http://mathworld.wolfram.com/CellularAutomaton.html">Mathworld,         Cellular Automaton</a></div>         <div>Rule 110 expands to 01110110</div>

    sage: test(interacts.geometry.unit_circle)
    Interactive function <function unit_circle at ...> with 2 widgets
      function: Dropdown(value=0, options=[('sin(x)', 0), ('cos(x)', 1), ('tan(x)', 2)], description=u'function')
      x: TransformFloatSlider(value=0.0, min=0.0, max=6.283185307179586, step=0.015707963267948967, description=u'x')
    <div style="white-space: normal;">Lines of the same color have         the same length</div>

    sage: test(interacts.geometry.trigonometric_properties_triangle)
    Interactive function <function trigonometric_properties_triangle at ...> with 3 widgets
      a0: IntSlider(value=30, min=0, max=360, step=1, description=u'A')
      a1: IntSlider(value=180, min=0, max=360, step=1, description=u'B')
      a2: IntSlider(value=300, min=0, max=360, step=1, description=u'C')
    <h2>Trigonometric Properties of a Triangle</h2>
    <script type="math/tex">\angle A = {60.0}^{\circ},</script> <script type="math/tex">\angle B = {45.0}^{\circ},</script> <script type="math/tex">\angle C = {75.0}^{\circ}</script>
    <script type="math/tex">AB = 1.93185165258,</script>  <script type="math/tex">BC = 1.73205080757,</script>  <script type="math/tex">CA = 1.41421356237</script>
    Area of triangle <script type="math/tex">ABC = 1.18301270189</script>

    sage: test(interacts.geometry.special_points)
    Interactive function <function special_points at ...> with 10 widgets
      title: HTMLText(value=u'<h2>Special points in triangle</h2>')
      a0: IntSlider(value=30, min=0, max=360, step=1, description=u'A')
      a1: IntSlider(value=180, min=0, max=360, step=1, description=u'B')
      a2: IntSlider(value=300, min=0, max=360, step=1, description=u'C')
      show_median: Checkbox(value=False, description=u'Medians')
      show_pb: Checkbox(value=False, description=u'Perpendicular Bisectors')
      show_alt: Checkbox(value=False, description=u'Altitudes')
      show_ab: Checkbox(value=False, description=u'Angle Bisectors')
      show_incircle: Checkbox(value=False, description=u'Incircle')
      show_euler: Checkbox(value=False, description=u"Euler's Line")

    sage: test(interacts.statistics.coin)
    Interactive function <function coin at ...> with 2 widgets
      n: IntSlider(value=1000, min=2, max=10000, step=100, description=u'Number of Tosses')
      interval: IntRangeSlider(value=(0, 0), min=0, max=1, step=1, description=u'Plotting range (y)')
    doctest:...: UserWarning: Attempting to set identical bottom==top results in singular transformations; automatically expanding.
    bottom=0.0, top=0.0
