from __future__ import absolute_import

from .api import Graphics3dRenderer

registered_renderers = {}

def register(renderer):
    global registered_renderers
    registered_renderers[renderer.__name__] = renderer

def unregister(renderer):
    del registered_renderers[renderer.__name__]
