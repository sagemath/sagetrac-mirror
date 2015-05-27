from __future__ import absolute_import

from .api import Graphics3dRenderer

registered_renderers = {}

def register(renderer):
    global registered_renderers
    registered_renderers[renderer.name] = renderer()

def unregister(renderer):
    global registered_renderers
    del registered_renderers[renderer.name]



# register the renderers in the sage library:
from .canvas3d import Canvas3dRenderer
from .jmol import JmolRenderer
from .obj import ObjRenderer
from .tachyon import TachyonRenderer
from .wavefront import WavefrontRenderer
from .x3d import X3dRenderer

register(Canvas3dRenderer)
register(JmolRenderer)
register(ObjRenderer)
register(TachyonRenderer)
register(WavefrontRenderer)
register(X3dRenderer)


"""
a test for java

        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        if not jdata.is_jvm_available():
            # We can only use JMol to generate preview if a jvm is installed
"""
