from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer


class ObjRenderer(Graphics3dRenderer):

    # grob = (gr)aphics (ob)ject
    def render_graphics3d(self, grob, render_params):
        """
        Unless otherwise changed, all rendering methods fall back to this
        one.
        """
        return ''

    # def render_graphics3d_group(self, grob, render_params):
    #     return self.render_graphics3d(grob, render_params)
    def render_transform_group(self, grob, render_params):
        return self.render_graphics3d_group(grob, render_params)

    def render_primitive_grobect(self, grob, render_params):
        return self.render_graphics3d(grob, render_params)
    def render_line(self, grob, render_params):
        return self.render_primitive_grobectd(grob, render_params)
    def render_point(self, grob, render_params):
        return self.render_primitive_grobect(grob, render_params)

    def render_index_face_set(self, grob, render_params):
        """
        Return an obj representation for ``grob``.

        TESTS::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Cylinder(1,1)
            sage: s = S.obj_repr(S.default_render_params())
        """
        transform = render_params.transform
        if transform is None:
            transform_point = tuple
        else:
            transform_point = transform.transform_point

        off = render_params.obj_vertex_offset

        grob_faces = grob.faces()

        points = ["v %g %g %g"%(transform_point(v)) for v in grob.vertices()]

        faces = [format_obj_face(indexed_face, off) for indexed_face in grob.index_faces()]
        if not grob.is_enclosed():
            back_faces = [format_obj_face(indexed_face, off, order=-1) for indexed_face in grob.index_faces()]
        else:
            back_faces = []

        render_params.obj_vertex_offset += len(grob.vertices())


        return ["g " + render_params.unique_name('obj'),
                "usemtl " + grob.texture.id,
                points,
                faces,
                back_faces]

    def render_box(self, grob, render_params):
        return self.render_index_face_set(grob, render_params)

    def render_parametric_surface(self, grob, render_params):
        grob.triangulate()
        return self.render_index_face_set(grob, render_params)
    def render_sphere(self, grob, render_params):
        return self.render_parametric_surface(grob, render_params)
    def render_cylinder(self, grob, render_params):
        return self.render_parametric_surface(grob, render_params)
    def render_torus(self, grob, render_params):
        return self.render_parametric_surface(grob, render_params)
    def render_cone(self, grob, render_params):
        return self.render_parametric_surface(grob, render_params)
    def render_mobius_strip(self, grob, render_params):
        return self.render_parametric_surface(grob, render_params)

    def render_implicit_surface(self, grob, render_params):
        grob.triangulate()
        return self.render_index_face_set(grob, render_params)

register(ObjRenderer)

def format_obj_face(indexed_face, off, order=1):
    return "f "+" ".join("%d"%(j+off) for j in indexed_face[::order])
