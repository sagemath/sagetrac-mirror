
class Graphics3dRenderer(object):

    def render_3d_group(self, obj):

    def render_sphere(self, obj):
        """
        Get a sphere, give back a string that will render the sphere with xxxx
        """
        return self.render_parametric_surface(obj)

    def render_parametric_surface(self, obj):
        obj.triangulate()
        return self.render_index_face_set(obj)
