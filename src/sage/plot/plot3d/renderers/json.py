from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer

class JSONRenderer(Graphics3dRenderer):

#    def rich_repr_graphics3d(self, obj, **kwds):
#        pass

#    def render_graphics3d(self, obj):
#        pass

    def render_index_face_set(self, obj, render_params):
        pass

    def render_implicit_surface(self, obj, render_params):
        obj.triangulate()

    def render_parametric_surface(self, obj, render_params):
        obj.triangulate()



register(JSONRenderer)
