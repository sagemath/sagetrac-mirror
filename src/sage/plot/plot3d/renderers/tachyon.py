from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer


class TachyonRenderer(Graphics3dRenderer):

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
        return self.render_parametric_surface(obj, render_params)
    def render_cylinder(self, obj, render_params):
        return self.render_parametric_surface(obj, render_params)
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






