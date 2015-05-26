from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer
from sage.repl.rich_output.buffer import OutputBuffer
from sage.repl.rich_output.output_graphics3d  import OutputSceneJmol

class JMOLRenderer(Graphics3dRenderer):
    name = 'Jmol'
    output_type = OutputSceneJmol
    def rich_repr(self, obj, **kwds):
        """
        Rich Representation as JMol scene

        INPUT:

        Optional keyword arguments are passed to JMol.

        OUTPUT:

        Instance of
        :class:`sage.repl.rich_output.output_graphics3d.OutputSceneJmol`.

        EXAMPLES::

            sage: sphere()._rich_repr_jmol()
            OutputSceneJmol container
        """
        from sage.misc.temporary_file import tmp_dir
        root_dir = os.path.abspath(tmp_dir())
        scene_zip     = os.path.join(root_dir, 'scene.spt.zip')
        preview_png   = os.path.join(root_dir, 'preview.png')
        opts = obj._process_viewing_options(kwds)
        zoom = opts['zoom']
        T = self._prepare_for_jmol(
            obj,
            opts['frame'],
            opts['axes'],
            opts['frame_aspect_ratio'],
            opts['aspect_ratio'],
            zoom,
        )
        self.export_jmol(T, scene_zip, **opts)
        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        # Java needs absolute paths
        # On cygwin, they should be native ones
        scene_native = scene_zip
        import sys
        if sys.platform == 'cygwin':
            from subprocess import check_output, STDOUT
            scene_native = check_output(['cygpath', '-w', scene_native],
                                        stderr=STDOUT).rstrip()
        script = '''set defaultdirectory "{0}"\nscript SCRIPT\n'''.format(scene_native)
        jdata.export_image(targetfile=preview_png, datafile=script,
                           image_type="PNG",
                           figsize=opts['figsize'][0])
        from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
        from sage.repl.rich_output.buffer import OutputBuffer
        scene_zip     = OutputBuffer.from_file(scene_zip)
        preview_png   = OutputBuffer.from_file(preview_png)
        return OutputSceneJmol(scene_zip, preview_png)

    def _prepare_for_jmol(self, obj, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        from sage.plot.plot import EMBEDDED_MODE
        if EMBEDDED_MODE:
            s = 6
        else:
            s = 3
        box_min, box_max = obj._rescale_for_frame_aspect_ratio_and_zoom(s, frame_aspect_ratio, zoom)
        a_min, a_max = obj._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return obj._transform_to_bounding_box(box_min, box_max, a_min, a_max, frame=frame,
                                               axes=axes, thickness=1,
                                               labels = True)   # jmol labels are implemented

    def export_jmol(self, obj, filename='jmol_shape.jmol', force_reload=False,
                    zoom=1, spin=False, background=(1,1,1), stereo=False,
                    mesh=False, dots=False,
                    perspective_depth = True,
                    orientation = (-764,-346,-545,76.39), **ignored_kwds):
                    # orientation chosen to look same as tachyon
        """
        A jmol scene consists of a script which refers to external files.
        Fortunately, we are able to put all of them in a single zip archive,
        which is the output of this call.

        EXAMPLES::

            sage: out_file = tmp_filename(ext=".jmol")
            sage: G = sphere((1, 2, 3), 5) + cube() + sage.plot.plot3d.shapes.Text("hi")
            sage: G.export_jmol(out_file)
            sage: import zipfile
            sage: z = zipfile.ZipFile(out_file)
            sage: z.namelist()
            ['obj_...pmesh', 'SCRIPT']

            sage: print z.read('SCRIPT')
            data "model list"
            2
            empty
            Xx 0 0 0
            Xx 5.5 5.5 5.5
            end "model list"; show data
            select *
            wireframe off; spacefill off
            set labelOffset 0 0
            background [255,255,255]
            spin OFF
            moveto 0 -764 -346 -545 76.39
            centerAt absolute {0 0 0}
            zoom 100
            frank OFF
            set perspectivedepth ON
            isosurface sphere_1  center {1.0 2.0 3.0} sphere 5.0
            color isosurface  [102,102,255]
            pmesh obj_... "obj_...pmesh"
            color pmesh  [102,102,255]
            select atomno = 1
            color atom  [102,102,255]
            label "hi"
            isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;

            sage: print z.read(z.namelist()[0])
            24
            0.5 0.5 0.5
            -0.5 0.5 0.5
            ...
            -0.5 -0.5 -0.5
            6
            5
            0
            1
            ...
        """
        from cStringIO import StringIO
        import zipfile
        from ..base import flatten_list
        render_params = obj.default_render_params()
        render_params.mesh = mesh
        render_params.dots = dots
        render_params.output_file = filename
        render_params.force_reload = render_params.randomize_counter = force_reload
        render_params.output_archive = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED, True)
        # Render the data
        all = flatten_list([obj.render(render_params,self), ""])

        f = StringIO()

        if render_params.atom_list:
            # Load the atom model
            f.write('data "model list"\n')
            f.write('%s\nempty\n' % (len(render_params.atom_list) + 1))
            for atom in render_params.atom_list:
                f.write('Xx %s %s %s\n' % atom)
            f.write('Xx 5.5 5.5 5.5\n') # so the zoom fits the box
            f.write('end "model list"; show data\n')
            f.write('select *\n')
            f.write('wireframe off; spacefill off\n')
            f.write('set labelOffset 0 0\n')


        # Set the scene background color
        f.write('background [%s,%s,%s]\n'%tuple([int(a*255) for a in background]))
        if spin:
            f.write('spin ON\n')
        else:
            f.write('spin OFF\n')
        if stereo:
            if stereo is True: stereo = "redblue"
            f.write('stereo %s\n' % stereo)
        if orientation:
            f.write('moveto 0 %s %s %s %s\n'%tuple(orientation))

        f.write('centerAt absolute {0 0 0}\n')
        f.write('zoom {0}\n'.format(zoom * 100))
        f.write('frank OFF\n') # jmol logo

        if perspective_depth:
            f.write('set perspectivedepth ON\n')
        else:
            f.write('set perspectivedepth OFF\n')

        # Put the rest of the object in
        f.write("\n".join(all))
        # Make sure the lighting is correct
        f.write("isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;\n")

        render_params.output_archive.writestr('SCRIPT', f.getvalue())
        render_params.output_archive.close()


    def render_graphics3d_image(self,obj):
        scene = obj._rich_repr_jmol(**opts)
        scene.preview_png.save_as(filename)

    def render_graphics3d(self, obj, render_params):
        return ''

    # def render_graphics3d_group(self, obj, render_params):
    #     return [g.render(render_params, renderer=self) for g in obj.all]
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
        Return a jmol representation for ``obj``.

        TESTS::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Cylinder(1,1)
            sage: S.show(viewer='jmol')   # indirect doctest
        """
        transform = render_params.transform
        vertices = obj.vertices()
        faces = obj.faces()

        obj._seperate_creases(render_params.crease_threshold)

        if transform is None:
            points = ["%g %g %g"%tuple(v) for v in vertices]
        else:
            points = ["%g %g %g"%transform.transform_point(v) for v in vertices]


        pmesh_faces = []
        for f in faces:
            if len(f) <= 4:
                s = "%d\n"%(len(f)+1)
                s += "\n".join("%d"%index for index in f.iter_index())
                s += "\n%d"%f.get_index(0)
                pmesh_faces.append(s)
            else:
                # stupid trianglulation
                f0 = f.get_index(0)
                for i in range(1, len(f)-1):
                    pmesh_faces.append("%d\n%d\n%d\n%d\n%d"%(4,f0,f.get_index(i),f.get_index(i+1),f0))

        # activation of coloring in jmol
        #if obj.global_texture:
        #    pmesh_faces = [format_pmesh_face(f, 1) for f in faces]
        #else:
        #    pmesh_faces = [format_pmesh_face(f, -1) for f in faces]

        # If a face has more than 4 vertices, it gets chopped up in
        # format_pmesh_face
        extra_faces = 0
        for f in faces:
            if len(f) >= 5:
                extra_faces += len(f)-3

        all = [str(len(vertices)),
               points,
               str(len(faces) + extra_faces),
               pmesh_faces]

        from ..base import flatten_list
        name = render_params.unique_name('obj')
        all = flatten_list(all)
        if render_params.output_archive:
            filename = "%s.pmesh" % (name)
            render_params.output_archive.writestr(filename, '\n'.join(all))
        else:
            filename = "%s-%s.pmesh" % (render_params.output_file, name)
            f = open(filename, 'w')
            for line in all:
                f.write(line)
                f.write('\n')
            f.close()

        if True:#obj.global_texture:
            s = 'pmesh {} "{}"\n{}'.format(name, filename,
                                           obj.texture.jmol_str("pmesh"))
        else:
            s = 'pmesh {} "{}"'.format(name, filename)

        # Turn on display of the mesh lines or dots?
        if render_params.mesh:
            s += '\npmesh %s mesh\n' % name
        if render_params.dots:
            s += '\npmesh %s dots\n' % name
        return s

    def render_box(self, obj, render_params):
        return self.render_index_face_set(obj, render_params)

    def render_parametric_surface(self, obj, render_params):
        obj.triangulate()
        return self.render_index_face_set(obj, render_params)
    def render_sphere(self, obj, render_params):
        r"""
        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import Sphere

        Jmol has native code for handling spheres::

            sage: S = Sphere(2)
            sage: S.jmol_repr(S.default_render_params())
            ['isosurface sphere_1  center {0 0 0} sphere 2.0\ncolor isosurface  [102,102,255]']
            sage: S.translate(10, 100, 1000).jmol_repr(S.default_render_params())
            [['isosurface sphere_1  center {10.0 100.0 1000.0} sphere 2.0\ncolor isosurface  [102,102,255]']]

        It cannot natively handle ellipsoids::

            sage: Sphere(1).scale(2, 3, 4).jmol_repr(S.testing_render_params())
            [['pmesh obj_2 "obj_2.pmesh"\ncolor pmesh  [102,102,255]']]

        Small spheres need extra hints to render well::

            sage: Sphere(.01).jmol_repr(S.default_render_params())
            ['isosurface sphere_1 resolution 100 center {0 0 0} sphere 0.01\ncolor isosurface  [102,102,255]']
        """
        from math import sqrt
        name = render_params.unique_name('sphere')
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
        if rad < 0.5:
            res = "resolution %s" % min(int(7/rad), 100)
        else:
            res = ""
        return "isosurface %s %s center {%s %s %s} sphere %s\n%s" % (name, res, cen[0], cen[1], cen[2], rad, obj.texture.jmol_str("isosurface"))


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




register(JMOLRenderer)
