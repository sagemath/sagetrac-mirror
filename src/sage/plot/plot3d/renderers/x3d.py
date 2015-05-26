from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer

class X3DRenderer(Graphics3dRenderer):
    def render_graphics3d_image(self,obj):
        scene = obj._rich_repr_jmol(**opts)
        scene.preview_png.save_as(filename)

    def render_index_face_set(self, obj, render_params):
        pass

    def render_implicit_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)

    def render_parametric_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)



register(X3DRenderer)
