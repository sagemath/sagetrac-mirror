r""" Enable/Disable automatic variable creation in sage.
This works by catching NameErrors
and defining the caught names as variables through `var`.

Beware that there is absolutely no lexical analysis of the code,
as differentiating between a typo and an intended variable is fairly hard.

Use wisely.

AUTHORS:

- Konstantin Kliakhandler (2016-11-24 sagedays79): Initial version
"""

import re
import sys,inspect


class __auto_vars(object):
  """ Enable/Disable automatic variable creation in sage.

  This works by catching NameErrors
  and defining the caught names as variables through `var`.

  Beware that there is absolutely no lexical analysis of the code,
  as differentiating between a typo and an intended variable is fairly hard.

  Use wisely.

  EXAMPLES::

    sage: auto_vars()         # not tested

    sage: a                   # not tested
    a

    sage: f = x^y - z/t       # not tested
    sage: f                   # not tested
    x^y - z/t

    sage: ## Note that typos are not detected, by design.
    sage: soolve()            # not tested
    soolve

    sage: def foo():          # not tested
    ....:    'Variables are also created if something is undefined inside a function!'
    ....:     return undefined
    ....:

    sage: foo()               # not tested
    undefined
  """
  def __init__(self):
    self.var_pattern = re.compile("'[a-zA-Z_][a-zA-Z_0-9]*'")
    self.enabled = False

  def __call__(self, state=True):
    """
    Enable/Disable automatic variable creation.
    `state`: Enabled/Disabled state.
    """
    self.enabled = state
    if state:
      ip = get_ipython()
      ip.set_custom_exc((NameError,), ipython_exception_wrapper)


auto_vars = __auto_vars()

def ipython_exception_wrapper(self, etype, value, tb, tb_offset=None):
  if not (etype is NameError and auto_vars.enabled):
      return tb

  m = value.args[0]
  result = auto_vars.var_pattern.search(m)
  if result is None:
    return tb

  name = result.group().strip("'")
  name_eval_src = "var('%s');" % name

  ip = get_ipython()
  try:
    ip.run_cell(name_eval_src, store_history=False, silent=True)

    ip.run_cell(inspect.getsource(sys.exc_info()[2].tb_next))
  except:
    ip.showtraceback()

  return None
