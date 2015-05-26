from __future__ import absolute_import

from .api import Graphics3dRenderer

registered_renderers = {}

def register(renderer):
    global registered_renderers
    registered_renderers[renderer.name] = renderer()

def unregister(renderer):
    del registered_renderers[renderer.name]




"""
a test for java

        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        if not jdata.is_jvm_available():
            # We can only use JMol to generate preview if a jvm is installed
"""
