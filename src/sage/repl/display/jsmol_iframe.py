"""
HTML Generator for JSmol

This is all an evil iframe hack to get JSmol to display 3-d graphics
while separating JSmol's j2s machinery from your actual web page.

There are some caveats for how to load JSmol, in particular it cannot
just load its code from a ``file://`` uri. To use a html file
generated by this module, you need

* A web server,

* The JSmol directory tree must be served by your web server,

* The output of :meth:`JSMolHtml.inner_html` or
  :meth:`JSMolHtml.outer_html` must be served by the same web server.

See https://github.com/phetsims/molecule-polarity/issues/6 for a
discussion of loading JSMol.
"""

import io
import os
import zipfile

from sage.cpython.string import bytes_to_str
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method


INNER_HTML_TEMPLATE = \
"""
<html>
<head>
  <style>
    * {{
      margin: 0;
      padding: 0;
      overflow: hidden;
    }}
    body, html {{      
      height: 100%;
      width: 100%;
    }}
  </style>
  <script type="text/javascript" src="{path_to_jsmol}/JSmol.min.js"></script>
</head>
<body>
  <script type="text/javascript">
    var script = {script};
    var Info = {{
      width: '{width}',
      height: '{height}',
      debug: false,
      disableInitialConsole: true,   // very slow when used with inline mesh
      color: '#3131ff',
      addSelectionOptions: false,
      use: 'HTML5',
      j2sPath: '{path_to_jsmol}/j2s',
      script: script,
    }};
    var jmolApplet0 = Jmol.getApplet('jmolApplet0', Info);
  </script>
</body>
</html>
"""


IFRAME_TEMPLATE = \
"""
<iframe srcdoc="{escaped_inner_html}" 
        width="{width}"
        height="{height}"
        style="border: 0;">
</iframe>
"""


OUTER_HTML_TEMPLATE = \
"""
<html>
<head>
  <title>JSmol 3D Scene</title>
</head>
</body>
  {iframe}
</body>
</html>
""".format(iframe=IFRAME_TEMPLATE)


class JSMolHtml(SageObject):

    def __init__(self, jmol, path_to_jsmol=None, width='100%', height='100%'):
        """
        INPUT:

        - ``jmol`` -- 3-d graphics or
          :class:`sage.repl.rich_output.output_graphics3d.OutputSceneJmol`
          instance. The 3-d scene to show.

        - ``path_to_jsmol`` -- string (optional, default is
          ``'/nbextensions/jsmol'``). The path (relative or absolute)
          where ``JSmol.min.js`` is served on the web server. 

        - ``width`` -- integer or string (optional, default:
          ``'100%'``). The width of the JSmol applet using CSS
          dimensions.

        - ``height`` -- integer or string (optional, default:
          ``'100%'``). The height of the JSmol applet using CSS
          dimensions.

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: JSMolHtml(sphere(), width=500, height=300)
            JSmol Window 500x300
        """
        from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
        if not isinstance(jmol, OutputSceneJmol):
            jmol = jmol._rich_repr_jmol()
        self._jmol = jmol
        self._zip = zipfile.ZipFile(io.BytesIO(self._jmol.scene_zip.get()))
        if path_to_jsmol is None:
            self._path = os.path.join('/', 'nbextensions', 'jsmol')
        else:
            self._path = path_to_jsmol
        self._width = width
        self._height = height

    @cached_method
    def script(self):
        r"""
        Return the JMol script file.

        This method extracts the Jmol script from the Jmol spt file (a
        zip archive) and inlines meshes.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: jsmol = JSMolHtml(OutputSceneJmol.example(), width=500, height=300)
            sage: jsmol.script()
            'data "model list"\n10\nempt...aliasdisplay on;\n'
        """
        script = []
        with self._zip.open('SCRIPT') as SCRIPT:
            for line in SCRIPT:
                if line.startswith(b'pmesh'):
                    command, obj, meshfile = line.split(b' ', 3)
                    assert command == b'pmesh'
                    if meshfile not in [b'dots\n', b'mesh\n']:
                        assert (meshfile.startswith(b'"') and
                                meshfile.endswith(b'"\n'))
                        meshfile = bytes_to_str(meshfile[1:-2])  # strip quotes
                        script += [
                            'pmesh {0} inline "'.format(bytes_to_str(obj)),
                            bytes_to_str(self._zip.open(meshfile).read()),
                            '"\n'
                        ]
                        continue
                script += [bytes_to_str(line)]
        return ''.join(script)

    def js_script(self):
        r"""
        The :meth:`script` as Javascript string.

        Since the many shortcomings of Javascript include multi-line
        strings, this actually returns Javascript code to reassemble
        the script from a list of strings.
        
        OUTPUT:

        String. Javascript code that evaluates to :meth:`script` as
        Javascript string.

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: jsmol = JSMolHtml(OutputSceneJmol.example(), width=500, height=300)
            sage: print(jsmol.js_script())
            [
              'data "model list"',
              ...
              'isosurface fullylit; pmesh o* fullylit; set antialiasdisplay on;',
            ].join('\n');
        """
        script = [r"["]
        for line in self.script().splitlines():
            script += [r"  '{0}',".format(line)]
        script += [r"].join('\n');"]
        return '\n'.join(script)
        
    def _repr_(self):
        """
        Return as string representation

        OUTPUT:
        
        String.

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: JSMolHtml(OutputSceneJmol.example(), width=500, height=300)._repr_()
            'JSmol Window 500x300'
        """
        return 'JSmol Window {0}x{1}'.format(self._width, self._height)

    def inner_html(self):
        """
        Return a HTML document containing a JSmol applet

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: jmol = JSMolHtml(OutputSceneJmol.example(), width=500, height=300)
            sage: print(jmol.inner_html())
            <html>
            <head>
              <style>
                * {
                  margin: 0;
                  padding: 0;
                    ...
            </html>
        """
        return INNER_HTML_TEMPLATE.format(
            script=self.js_script(),
            width=self._width,
            height=self._height,
            path_to_jsmol=self._path,
        )
        
    def iframe(self):
        """
        Return HTML iframe

        OUTPUT:
        
        String.
        
        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: jmol = JSMolHtml(OutputSceneJmol.example())
            sage: print(jmol.iframe())
            <iframe srcdoc="
            ...
            </iframe>
        """
        escaped_inner_html = self.inner_html().replace('"', '&quot;')
        iframe = IFRAME_TEMPLATE.format(
            script=self.js_script(),
            width=self._width,
            height=self._height,
            escaped_inner_html=escaped_inner_html,
        )
        return iframe

    def outer_html(self):
        """
        Return a HTML document containing an iframe with a JSmol applet

        OUTPUT:

        String

        EXAMPLES::

            sage: from sage.repl.display.jsmol_iframe import JSMolHtml
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: jmol = JSMolHtml(OutputSceneJmol.example(), width=500, height=300)
            sage: print(jmol.outer_html())
            <html>
            <head>
              <title>JSmol 3D Scene</title>
            </head>
            </body>
            <BLANKLINE>
            <iframe srcdoc="
                    ...
            </html>
        """
        escaped_inner_html = self.inner_html().replace('"', '&quot;')
        outer = OUTER_HTML_TEMPLATE.format(
            script=self.js_script(),
            width=self._width,
            height=self._height,
            escaped_inner_html=escaped_inner_html,
        )
        return outer
