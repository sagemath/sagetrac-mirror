from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer

from sage.repl.rich_output.output_catalog import OutputSceneCanvas3d

class Canvas3dRenderer(Graphics3dRenderer):
    name='canvas3d'
    output_types=(OutputSceneCanvas3d,)
    def rich_repr(self, obj, **kwds):
        r"""
        Rich Representation as Canvas3D Scene

        INPUT:

        Optional keyword arguments.

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneCanvas3d`.

        EXAMPLES::

            sage: out = sphere()._rich_repr_canvas3d()
            sage: out
            OutputSceneCanvas3d container
            sage: out.canvas3d.get()
            "[{vertices:[{x:0,y:0,z:-1},...,color:'#6666ff'}]"
        """
        opts = obj._process_viewing_options(kwds)
        aspect_ratio = opts['aspect_ratio'] # this necessarily has a value now
        frame_aspect_ratio = opts['frame_aspect_ratio']
        zoom = opts['zoom']
        frame = opts['frame']
        axes = opts['axes']
        T = obj._prepare_for_render(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
        data = flatten_list(T.render(T.default_render_params(), renderer=self))
        canvas3d = '[' + ','.join(data) + ']'
        return OutputSceneCanvas3d(canvas3d)
        
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
        Return a canvas3d representation for ``self``.

        TESTS:

        A basic test with a triangle::

            sage: G = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: G.json_repr(G.default_render_params())
            ["{vertices:[{x:0,y:0,z:1},{x:1,y:1,z:1},{x:2,y:0,z:1}],faces:[[0,1,2]],color:'#0000ff'}"]

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.json_repr(S.default_render_params())
            ["{vertices:[{x:2,y:0,z:0},{x:0,y:2,z:0},{x:0,y:0,z:2},{x:0,y:1,z:1},{x:1,y:0,z:1},{x:1,y:1,z:0}],faces:[[0,4,5],[3,4,5],[2,3,4],[1,3,5]],face_colors:['#ff0000','#ff9900','#cbff00','#33ff00']}"]
        """
        transform = render_params.transform
        vertices = obj.vertices()
        faces = obj.faces()

        if transform is None:
            vertices_str = "[{}]".format(",".join(["{x:%g,y:%g,z:%g}"%tuple(v) for v in vertices]))
        else:
            vertices_str = "[{}]".format(",".join(
                ["{x:%g,y:%g,z:%g}"%transform.transform_point(v) for v in vertices]
            ))

        faces_str = "[{}]".format(",".join(
            ["[{}]".format(",".join(
                [str(v) for v in face.iter_index()]
            )) for face in faces]
        ))

        if True:#obj.global_texture:
            color_str = "'#{}'".format(obj.texture.hex_rgb())
            return "{vertices:%s,faces:%s,color:%s}"%(vertices_str, faces_str, color_str)
        # else:
        #     color_str = "[{}]".format(",".join(["'{}'".format(
        #             Color(obj._faces[i].color.r,
        #                   obj._faces[i].color.g,
        #                   obj._faces[i].color.b).html_color())
        #                                     for i from 0 <= i < obj.fcount]))
        #     return "{vertices:%s,faces:%s,face_colors:%s}"%(vertices_str, faces_str, color_str)


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

register(Canvas3dRenderer)


# cdef inline format_json_vertex(point_c P):
#     cdef char ss[100]
#     cdef Py_ssize_t r = sprintf_3d(ss, "{x:%g,y:%g,z:%g}", P.x, P.y, P.z)
#     return PyString_FromStringAndSize(ss, r)

# cdef inline format_json_face(face_c face):
#     return "[{}]".format(",".join([str(face.vertices[i])
#                                    for i from 0 <= i < face.n]))

