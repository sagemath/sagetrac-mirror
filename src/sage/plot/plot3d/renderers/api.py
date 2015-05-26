
class Graphics3dRenderer(object):

    def render_graphics3d(self, obj):
        raise NotImplementedError

    def render_graphics3d_group(self, obj):
        return self.render_graphics3d(obj)
    def render_transform_group(self, obj):
        return self.render_graphics3d_group(obj)

    def render_primitive_object(self, obj):
        return self.render_graphics3d(obj)
    # def render_line(self, obj):
    #     return self.render_primitive_objectd(obj)
    # def render_point(self, obj):
    #     return self.render_primitive_object(obj)

    def render_index_face_set(self, obj):
        return self.render_graphics3d(obj)
    # def render_box(self, obj):
    #     return self.render_index_face_set(obj)

    def render_parametric_surface(self, obj):
        obj.triangulate()
        return self.render_index_face_set(obj)
    # def render_sphere(self, obj):
    #     return self.render_parametric_surface(obj)
    # def render_cylinder(self, obj):
    #     return self.render_parametric_surface(obj)
    # def render_torus(self, obj):
    #     return self.render_parametric_surface(obj)
    # def render_cone(self, obj):
    #     return self.render_parametric_surface(obj)
    # def render_mobius_strip(self, obj):
    #     return self.render_parametric_surface(obj)

    def render_implicit_surface(self, obj):
        obj.triangulate()
        return self.render_index_face_set(obj)

    def render_triangle_plot(self, obj): 
        obj.triangulate()
        return self.render_index_face_set(obj) #depends on TrianglePlot as subclass of IFS

