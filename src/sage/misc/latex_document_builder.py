from sage.misc.temporary_file import tmp_dir, graphics_filename
from sage.misc.viewer import pdf_viewer
import subprocess

COMMON_HEADER = \
r'''\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{amsfonts}
\usepackage{graphicx}
\usepackage{mathrsfs}
\pagestyle{empty}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
'''

LATEX_HEADER = (
r'''\documentclass{article}

''' + COMMON_HEADER +
r'''
\oddsidemargin 0.0in
\evensidemargin 0.0in
\textwidth 6.45in
\topmargin 0.0in
\headheight 0.0in
\headsep 0.0in
\textheight 9.0in
''')

class LatexDocument:
    class sage_object:
        def __init__(self, obj):
            self.obj = obj
        def getTex( self ):
            res = "\n\\["
            res += ( "\n" + latex( self.obj ) )
            res += "\n\\]"
            return res

    class latex:
        def __init__( self, tex ):
            self.tex = tex
        def getTex( self ):
            return self.tex

    def __init__(self, title = "Sage"):
        self.document = []
        self.set_title( title )

    def get_tex( self ):
        res = ""
        res += self.get_preambule()
        res += "\n\\begin{document}\n"
        res += "\n\\begin{center}{\\Large\\bf %s}\\end{center}\n"%(self.title)
        res += self.get_document()
        res += "\n\n\\end{document}"
        return res

    def export_tex( self, file_path ):
        open(file_path,'w').write( self.get_tex() )

    def export_pdf( self, file_path ):
        tmp = tmp_dir()
        tex_file = os.path.join(tmp, "sage.tex")
        self.export_tex( tex_file )
        subprocess.call(
            ["pdflatex", tex_file] 
            , stdout=subprocess.PIPE, stderr=subprocess.PIPE
            , cwd=tmp
        )
        pdf_file = os.path.join(tmp, "sage.pdf")
        subprocess.call(
            ["cp", pdf_file, file_path]
            , stdout=subprocess.PIPE, stderr=subprocess.PIPE
            , cwd=tmp
        )

    def set_title( self, title ):
        self.title = title

    def get_preambule( self ):
        res = ""
        res += ('\n' + LATEX_HEADER)
        return res

    def get_document( self ):
        res = ""
        for i in self.document:
            res += i.getTex()
        return res

    def add_latex( self, tex ):
        self.document.append( self.latex( tex ) )

    def add( self, sage_object ):
        self.document.append( self.sage_object( sage_object ) )

    def __repr__( self ):
        return self.get_document()

    def add_package_to_preamble_if_available(self, package_name):
        assert isinstance(package_name, str)
        if self.has_file(package_name+".sty"):
            self.add_to_preamble("\\usepackage{%s}\n"%package_name)

    def view( self ):
        tmp = tmp_dir('sage_viewer')
        pdf_file = os.path.join(tmp, "sage.pdf")
        self.export_pdf( pdf_file )
        subprocess.call(
            '%s %s'%(pdf_viewer(), pdf_file), 
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

