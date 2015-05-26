from __future__ import absolute_import

import os
from . import register, Graphics3dRenderer

class JMOLRenderer(Graphics3dRenderer):

    def rich_repr_graphics3d(self, obj, **kwds):
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
        T = obj._prepare_for_jmol(
            opts['frame'],
            opts['axes'],
            opts['frame_aspect_ratio'],
            opts['aspect_ratio'],
            zoom,
        )
        T.export_jmol(scene_zip, **opts)
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


    def render_graphics3d_image(self,obj):
        scene = obj._rich_repr_jmol(**opts)
        scene.preview_png.save_as(filename)

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

        sig_on()
        if transform is None:
            points = ["%g %g %g"%v for v in vertices]
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

        from base import flatten_list
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

        if obj.global_texture:
            s = 'pmesh {} "{}"\n{}'.format(name, filename,
                                           obj.texture.jmol_str("pmesh"))
        else:
            s = 'pmesh {} "{}"'.format(name, filename)

        # Turn on display of the mesh lines or dots?
        if render_params.mesh:
            s += '\npmesh %s mesh\n' % name
        if render_params.dots:
            s += '\npmesh %s dots\n' % name
        return [s]



    def render_implicit_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)

    def render_parametric_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)





register(JMOLRenderer)
