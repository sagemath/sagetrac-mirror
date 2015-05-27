from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer


# Needs to go somewhere:
"""
The tachyon representation of a group is simply the concatenation of
the representations of its objects.

EXAMPLES::

    sage: G = sphere() + sphere((1,2,3))
    sage: G.tachyon_repr(G.default_render_params())
    [['Sphere center 0.0 0.0 0.0 Rad 1.0 texture...'],
     ['Sphere center 1.0 2.0 3.0 Rad 1.0 texture...']]
"""
"""
Transformations for Tachyon are applied at the leaf nodes.

EXAMPLES::

    sage: G = sphere((1,2,3)).scale(2)
    sage: G.tachyon_repr(G.default_render_params())
    [['Sphere center 2.0 4.0 6.0 Rad 2.0 texture...']]
"""


class TachyonRenderer(Graphics3dRenderer):
    name = 'Tachyon'
    def render_graphics3d(self, obj, render_params):
        """
        Unless otherwise changed, all rendering methods fall back to this
        one.
        """
        return ''

    # def render_graphics3d_group(self, obj, render_params):
    #     return self.render_graphics3d(obj, render_params)
    # def render_transform_group(self, obj, render_params):
    #     return self.render_graphics3d_group(obj, render_params)

    def render_primitive_object(self, obj, render_params):
        return self.render_graphics3d(obj, render_params)
    def render_line(self, obj, render_params):
        return self.render_primitive_objectd(obj, render_params)
    def render_point(self, obj, render_params):
        return self.render_primitive_object(obj, render_params)

    def render_index_face_set(self, obj, render_params):

        """
        Return a tachyon object for ``self``.

        EXAMPLES:

        A basic test with a triangle::

            sage: G = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: s = G.tachyon_repr(G.default_render_params()); s
            ['TRI V0 0 0 1 V1 1 1 1 V2 2 0 1', ...]

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.tachyon_repr(S.default_render_params())
            ['TRI V0 2 0 0 V1 1 0 1 V2 1 1 0',
            'TEXTURE... AMBIENT 0.3 DIFFUSE 0.7 SPECULAR 0 OPACITY 1.0... COLOR 1 0 0 ... TEXFUNC 0',...]
        """
        transform = render_params.transform
        lines = []
        faces = obj.faces()

        if transform is None:
            transform_point = tuple
        else:
            transform_point = transform.transform_point

        for face in faces:
            P = transform_point(face[0])
            Q = transform_point(face[1])
            R = transform_point(face[2])
            lines.append(format_tachyon_triangle(P, Q, R))
            if True:#obj.global_texture:
                lines.append(obj.texture.id)
            else:
                lines.append(format_tachyon_texture(face.color))
            # stupid triangulation of polygon
            if len(face) > 3:
                for k in range(3,len(face)):
                    Q = R
                    R = transform_point(face[k])
                    lines.append(format_tachyon_triangle(P, Q, R))
                    if True:#obj.global_texture:
                        lines.append(obj.texture.id)
                    else:
                        lines.append(format_tachyon_texture(face.color))
        return lines


    def render_box(self, obj, render_params):
        return self.render_index_face_set(obj, render_params)

    def render_parametric_surface(self, obj, render_params):
        obj.triangulate()
        return self.render_index_face_set(obj, render_params)
    def render_sphere(self, obj, render_params):
        r"""
        Tachyon can natively handle spheres. Ellipsoids rendering is done
        as a parametric surface.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere
            sage: S = Sphere(2)
            sage: S.tachyon_repr(S.default_render_params())
            'Sphere center 0 0 0 Rad 2.0 texture...'
            sage: S.translate(1, 2, 3).scale(3).tachyon_repr(S.default_render_params())
            [['Sphere center 3.0 6.0 9.0 Rad 6.0 texture...']]
            sage: S.scale(1,1/2,1/4).tachyon_repr(S.default_render_params())
            [['TRI V0 0 0 -0.5 V1 0.308116 0.0271646 -0.493844 V2 0.312869 0 -0.493844',
              'texture...',
               ...
              'TRI V0 0.308116 -0.0271646 0.493844 V1 0.312869 0 0.493844 V2 0 0 0.5',
              'texture...']]
        """
        from math import sqrt
        transform = render_params.transform
        if not (transform is None or transform.is_uniform()):
            return self.render_parametric_surface(obj, render_params)

        if transform is None:
            cen = (0,0,0)
            rad = obj.radius
        else:
            cen = transform.transform_point((0,0,0))
            radv = transform.transform_vector((obj.radius,0,0))
            rad = sqrt(sum([x*x for x in radv]))
        return "Sphere center %s %s %s Rad %s %s" % (cen[0], cen[1], cen[2], rad, obj.texture.id)

    def render_cylinder(self, obj, render_params):
        r"""
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Cylinder
            sage: C = Cylinder(1/2, 4, closed=False)
            sage: C.tachyon_repr(C.default_render_params())
            'FCylinder\n   Base 0 0 0\n   Apex 0 0 4.0\n   Rad 0.5\n   texture...     '
            sage: C = Cylinder(1, 2)
            sage: C.tachyon_repr(C.default_render_params())
                ['Ring Center 0 0 0 Normal 0 0 1 Inner 0 Outer 1.0 texture...',
                 'FCylinder\n   Base 0 0 0\n   Apex 0 0 2.0\n   Rad 1.0\n   texture...     ',
                 'Ring Center 0 0 2.0 Normal 0 0 1 Inner 0 Outer 1.0 texture...']
        """
        transform = render_params.transform
        if not (transform is None or transform.is_uniform_on([(1,0,0),(0,1,0)])):
            # Tachyon can't do squashed
            return ParametricSurface.tachyon_repr(obj, render_params)

        base, top = obj.get_endpoints(transform)
        rad = obj.get_radius(transform)
        cyl = """FCylinder
   Base %s %s %s
   Apex %s %s %s
   Rad %s
   %s     """%(base[0], base[1], base[2], top[0], top[1], top[2], rad, obj.texture.id)
        if obj.closed:
            normal = (0,0,1)
            if transform is not None:
                normal = transform.transform_vector(normal)
            base_cap = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % (base[0], base[1], base[2], normal[0], normal[1], normal[2], rad, obj.texture.id)
            top_cap  = """Ring Center %s %s %s Normal %s %s %s Inner 0 Outer %s %s"""  \
                       % ( top[0],  top[1],  top[2], normal[0], normal[1], normal[2], rad, obj.texture.id)
            return [base_cap, cyl, top_cap]
        else:
            return cyl

    def render_torus(self, obj, render_params):
        return self.render_parametric_surface(obj, render_params)
    def render_cone(self, obj, render_params):
        return self.render_parametric_surface(obj, render_params)
    def render_mobius_strip(self, obj, render_params):
        return self.render_parametric_surface(obj, render_params)

    def render_implicit_surface(self, obj, render_params):
        obj.triangulate()
        return self.render_index_face_set(obj, render_params)

register(TachyonRenderer)


def format_tachyon_texture(rgb):
    return "TEXTURE\n AMBIENT 0.3 DIFFUSE 0.7 SPECULAR 0 OPACITY 1.0\n COLOR %g %g %g \n TEXFUNC 0\n" % (rgb.r, rgb.g, rgb.b)

def format_tachyon_triangle(P, Q, R):
    return "TRI V0 %g %g %g "%P + "V1 %g %g %g "%Q + "V2 %g %g %g\n"%R






