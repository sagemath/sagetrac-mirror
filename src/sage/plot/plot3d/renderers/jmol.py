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

    def repr_index_face_set(self, obj, render_params):
        pass

    def repr_implicit_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)

    def repr_parametric_surface(self, obj, render_params):
        obj.triangulate()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet.jmol_repr(obj, render_params)



register(JMOLRenderer)
