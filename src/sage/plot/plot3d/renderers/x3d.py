from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer


class X3dRenderer(Graphics3dRenderer):

    def render_graphics3d(self, obj, render_params):
        """
        Unless otherwise changed, all rendering methods fall back to this
        one.
        """
        return ''

    def render_graphics3d_group(self, obj, render_params):
        """
        The x3d representation of a group is simply the concatenation of
        the representation of its objects.

        EXAMPLES::

            sage: G = sphere() + sphere((1,2,3))
            sage: print G.x3d_str()
            <Transform translation='0 0 0'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
            <Transform translation='1 2 3'>
            <Shape><Sphere radius='1.0'/><Appearance><Material diffuseColor='0.4 0.4 1.0' shininess='1.0' specularColor='0.0 0.0 0.0'/></Appearance></Shape>
            </Transform>
        """
        return "\n".join([g.x3d_str() for g in obj.all])

    def render_transform_group(self, obj, render_params):
        return self.render_graphics3d_group(obj, render_params)

    def render_primitive_object(self, obj, render_params):
        return self.render_graphics3d(obj, render_params)
    def render_line(self, obj, render_params):
        return self.render_primitive_objectd(obj, render_params)
    def render_point(self, obj, render_params):
        return self.render_primitive_object(obj, render_params)

    def render_index_face_set(self, obj, render_params):
        return self.render_graphics3d(obj, render_params)
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

register(X3dRenderer)

