"""
Doc.

* Each subclass of Graphics3d has a render_method naming the method which should be used to render it.
* To implement a new renderer for 3d graphics, make a subclass of Graphics3dRenderer and implement rendering methods.
* To implement a new 3d graphics class, update the API here with a rendering method and fallback.  Update existing renderer methods as appropriate.

"""
class Graphics3dRenderer(object):

    def render_graphics3d(self, obj, render_params):
        """
        Unless otherwise changed, all rendering methods fall back to this
        one.
        """
        return ''

    def render_graphics3d_group(self, obj, render_params):
        """
        By default, apply render to each object in the group
        """
        return [g.render(render_params, renderer=self) for g in obj.all]
    def render_transform_group(self, obj, render_params):
        """
        Apply transformation to leaf nodes.
        """
        render_params.push_transform(obj.get_transformation())
        rep = [g.render(render_params, renderer=self) for g in obj.all]
        render_params.pop_transform()
        return rep

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


