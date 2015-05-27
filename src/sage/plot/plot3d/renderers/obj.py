from __future__ import absolute_import

import os
from .api import Graphics3dRenderer

# needs to go somewhere:
"""
The obj representation of a group is simply the concatenation of
the representation of its objects.

EXAMPLES::

    sage: G = tetrahedron() + tetrahedron().translate(10, 10, 10)
    sage: G.obj_repr(G.default_render_params())
    [['g obj_1',
      'usemtl ...',
      ['v 0 0 1',
       'v 0.942809 0 -0.333333',
       'v -0.471405 0.816497 -0.333333',
       'v -0.471405 -0.816497 -0.333333'],
      ['f 1 2 3', 'f 2 4 3', 'f 1 3 4', 'f 1 4 2'],
      []],
     [['g obj_2',
       'usemtl ...',
       ['v 10 10 11',
        'v 10.9428 10 9.66667',
        'v 9.5286 10.8165 9.66667',
        'v 9.5286 9.1835 9.66667'],
       ['f 5 6 7', 'f 6 8 7', 'f 5 7 8', 'f 5 8 6'],
       []]]]
"""
"""
Transformations for .obj files are applied at the leaf nodes.

EXAMPLES::

    sage: G = cube().scale(4).translate(1, 2, 3)
    sage: G.obj_repr(G.default_render_params())
    [[['g obj_1',
       'usemtl ...',
       ['v 3 4 5',
        'v -1 4 5',
        'v -1 0 5',
        'v 3 0 5',
        'v 3 4 1',
        'v -1 4 1',
        'v 3 0 1',
        'v -1 0 1'],
       ['f 1 2 3 4',
        'f 1 5 6 2',
        'f 1 4 7 5',
        'f 6 5 7 8',
        'f 7 4 3 8',
        'f 3 2 6 8'],
       []]]]
"""
class ObjRenderer(Graphics3dRenderer):
    name = 'Obj'
    # grob = (gr)aphics (ob)ject
    def render_graphics3d(self, grob, render_params):
        """
        Unless otherwise changed, all rendering methods fall back to this
        one.
        """
        return ''

    # def render_graphics3d_group(self, grob, render_params):
    #     return self.render_graphics3d(grob, render_params)
    # def render_transform_group(self, grob, render_params):
    #     return self.render_graphics3d_group(grob, render_params)

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
