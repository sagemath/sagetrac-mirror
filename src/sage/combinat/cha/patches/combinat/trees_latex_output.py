"""
Patches pour les sorties latex des arbres

TRAC : 14498

"""
from sage.combinat.cha.patches.monkey_patching import MonkeyPatch
from sage.combinat.abstract_tree import AbstractTree
from sage.misc.latex import latex
latex.add_package_to_preamble_if_available("tikz")
latex.add_to_mathjax_avoid_list("tikz")

class _(MonkeyPatch, AbstractTree):

    def _latex_(self):
        '''
        Nice output which can be easily modified

        TESTS::

            sage: latex(BinaryTree([[[],[]],[[],None]]))
            { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
            ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$$}
            ;}\newcommand{\nodec}{\node[draw,circle] (c) {$$}
            ;}\newcommand{\noded}{\node[draw,circle] (d) {$$}
            ;}\newcommand{\nodee}{\node[draw,circle] (e) {$$}
            ;}\newcommand{\nodef}{\node[draw,circle] (f) {$$}
            ;}\begin{tikzpicture}[auto]
            \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \&         \&         \& \nodea  \&         \&         \&         \\ 
                     \& \nodeb  \&         \&         \&         \& \nodee  \&         \\ 
             \nodec  \&         \& \noded  \&         \& \nodef  \&         \&         \\
            };
            <BLANKLINE>
            \path[ultra thick, red] (b) edge (c) edge (d)
                (e) edge (f)
                (a) edge (b) edge (e);
            \end{tikzpicture}}
        '''
        # latex environnement : TikZ
        begin_env = "\\begin{tikzpicture}[auto]\n"
        end_env = "\\end{tikzpicture}"
        # it uses matrix trick to place each node
        matrix_begin = "\\matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{\n"
        matrix_end = "\\\\\n};\n"
        # a basic path to each edges
        path_begin = "\\path[ultra thick, red] "
        path_end = ";\n"
        # to make a pretty output, it creates one LaTeX command for
        # each node
        cmd = "\\node"
        new_cmd1 = "\\newcommand{" + cmd
        new_cmd2 = "}{\\node[draw,circle] ("
        new_cmd3 = ") {$"
        new_cmd4 = "$}\n;}"
        # some variables to simplify code
        sep = "\\&"
        space = " " * 9
        sepspace = sep + space
        spacesep = space + sep
        node_to_str = lambda node: " " + node + " " * (
                    len(space) - 1 - len(node)
        )
        # # TODO:: modify how to create nodes --> new_cmd : \\node[...] in
        # # create_node
        num = [0]

        def resolve(self):
            nodes = []
            matrix = []
            edges = []

            def create_node(self):
                '''
                create a name (infixe reading)
                 -> ex: b
                create a new command:
                 -> ex: \newcommand{\nodeb}{\node[draw,circle] (b) {$$};
                return the name and the command to build:
                  . the matrix
                  . and the edges
                '''
                name = reduce(
                    lambda x, y: x + y,
                    map(
                        lambda x: chr(ord(x) + 49),
                        list(str(num[0]))),
                    "")
                node = cmd + name
                nodes.append((name,
                    (str(self.label()) if hasattr(self, "label") else ""))
                )
                num[0] += 1
                return node, name

            def empty_tree():
                '''
                TESTS::

                    sage: t = BinaryTree()
                    sage: print latex(t)
                    { \begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \\
                    };
                    \end{tikzpicture}}
                '''
                matrix.append(space)

            def one_node_tree(self):
                '''
                TESTS::

                    sage: t = BinaryTree([]); print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                    };
                    \end{tikzpicture}}
                    sage: t = OrderedTree([]); print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                    };
                    \end{tikzpicture}}
                '''
                node, _ = create_node(self)
                matrix.append(node_to_str(node))

            def concat_matrix(mat, mat2):
                lmat = len(mat)
                lmat2 = len(mat2)
                for i in range(max(lmat, lmat2)):
                    # mat[i] --> n & n & ...
                    # mat2[i] -> n' & n' & ...
                    # ==> n & n & ... & n' & n' & ...
                    try:
                        mat[i] += sep + mat2[i]
                    except:
                        if i >= lmat:
                            if i != 0:
                                # mat[i] doesn't exist but
                                # mat[0] has k "&"
                                # mat2[i] -> n' & n' & ...
                                # ==> (_ &)*k+1 n' & n' & ...
                                nb_of_and = mat[0].count(sep) - mat2[0].count(sep)
                                mat.append(spacesep * (nb_of_and) + mat2[i])
                            else:
                                # mat is empty
                                # mat2[i] -> n' & n' & ...
                                # ==> mat2
                                mat.extend(mat2)
                                return
                        else:
                            # mat[i] -> n & n & ...
                            # mat2[i] doesn't exist but mat2[0] exists
                            # # and has k "&"
                            # NOTE:: i != 0 because that is a no-empty subtree.
                            # ==> n & n & ... (& _)*k+1
                            nb_of_and = mat2[0].count(sep)
                            mat[i] += sepspace * (nb_of_and + 1)

            def tmp(subtree, edge, nodes, edges, matrix):
                if not subtree.is_empty():
                    # # create representation of the subtree
                    nodes_st, matrix_st, edges_st = resolve(subtree)
                    # # add its nodes to the "global" nodes set
                    nodes.extend(nodes_st)
                    # # create a new edge between the root and the subtree
                    edge.append(nodes_st[0][0])
                    # # add the subtree edges to the "global" edges set
                    edges.extend(edges_st)
                    # # build a new matrix by concatenation
                    concat_matrix(matrix, matrix_st)
                else:
                    concat_matrix(matrix, [space])

            def pair_nodes_tree(self, nodes, edges, matrix):
                """
                TESTS::

                    sage: t = OrderedTree([[[],[]],[[],[]]]).\
                    ....:    canonical_labelling(); print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\newcommand{\nodef}{\node[draw,circle] (f) {$6$}
                    ;}\newcommand{\nodeg}{\node[draw,circle] (g) {$7$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \&         \&         \& \nodea  \&         \&         \&         \\
                             \& \nodeb  \&         \&         \&         \& \nodee  \&         \\
                     \nodec  \&         \& \noded  \&         \& \nodef  \&         \& \nodeg  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d)
                        (e) edge (f) edge (g)
                        (a) edge (b) edge (e);
                    \end{tikzpicture}}
                    sage: t = BinaryTree([[],[[],[]]]); print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \&         \&         \\
                     \nodeb  \&         \&         \& \nodec  \&         \\
                             \&         \& \noded  \&         \& \nodee  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (c) edge (d) edge (e)
                        (a) edge (b) edge (c);
                    \end{tikzpicture}}
                """
                # build all subtree matrices.
                node, name = create_node(self)
                edge = [name]
                split = int(len(self) / 2)
                # the left part
                for i in range(split):
                    tmp(self[i], edge, nodes, edges, matrix)
                # # prepare the root line
                nb_of_and = matrix[0].count(sep)
                # the middle
                for i in range(len(matrix)):
                    matrix[i] += sepspace
                # the right part
                for i in range(split, len(self)):
                    tmp(self[i], edge, nodes, edges, matrix)

                # # create the root line
                root_line = (spacesep * (nb_of_and + 1) + node_to_str(node) +
                    sepspace * (matrix[0].count(sep) - nb_of_and - 1))
                matrix.insert(0, root_line)
                # add edges from the root
                edges.append(edge)

            def odd_nodes_tree(self, nodes, edges, matrix):
                '''
                TESTS::

                    sage: t = OrderedTree([[]]).canonical_labelling()
                    sage: print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                     \nodea  \\
                     \nodeb  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[]]]).canonical_labelling()
                    sage: print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \\
                             \& \nodeb  \&         \\
                     \nodec  \&         \& \noded  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d)
                        (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[],[]]]).canonical_labelling()
                    sage: print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \\
                             \& \nodeb  \&         \\
                     \nodec  \& \noded  \& \nodee  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d) edge (e)
                        (a) edge (b);
                    \end{tikzpicture}}
                    sage: t = OrderedTree([[[],[],[]],[],[]]).\
                    ....:     canonical_labelling(); print latex(t)
                    { \newcommand{\nodea}{\node[draw,circle] (a) {$1$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$2$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$3$}
                    ;}\newcommand{\noded}{\node[draw,circle] (d) {$4$}
                    ;}\newcommand{\nodee}{\node[draw,circle] (e) {$5$}
                    ;}\newcommand{\nodef}{\node[draw,circle] (f) {$6$}
                    ;}\newcommand{\nodeg}{\node[draw,circle] (g) {$7$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \&         \&         \& \nodea  \&         \\
                             \& \nodeb  \&         \& \nodef  \& \nodeg  \\
                     \nodec  \& \noded  \& \nodee  \&         \&         \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (b) edge (c) edge (d) edge (e)
                        (a) edge (b) edge (f) edge (g);
                    \end{tikzpicture}}
                '''
                # build all subtree matrices.
                node, name = create_node(self)
                edge = [name]
                split = int(len(self) / 2)
                # the left part
                for i in range(split):
                    tmp(self[i], edge, nodes, edges, matrix)
                # # prepare the root line
                if len(matrix) != 0:
                    nb_of_and = matrix[0].count(sep)
                    sizetmp = len(matrix[0])
                else:
                    nb_of_and = 0
                    sizetmp = 0
                # the middle
                tmp(self[split], edge, nodes, edges, matrix)
                nb_of_and += matrix[0][sizetmp:].split("node")[0].count(sep)

                # the right part
                for i in range(split + 1, len(self)):
                    tmp(self[i], edge, nodes, edges, matrix)

                # # create the root line
                root_line = (spacesep * (nb_of_and) + node_to_str(node) +
                    sepspace * (matrix[0].count(sep) - nb_of_and))
                matrix.insert(0, root_line)
                # add edges from the root
                edges.append(edge)
            if self.is_empty():
                empty_tree()
            elif len(self) == 0 or all(subtree.is_empty() for subtree in self):
                one_node_tree(self)
            elif len(self) % 2 == 0:
                pair_nodes_tree(self, nodes, edges, matrix)
            else:
                odd_nodes_tree(self, nodes, edges, matrix)
            return nodes, matrix, edges

        nodes, matrix, edges = resolve(self)

        def make_cmd(nodes):
            cmds = []
            for name, label in nodes:
                cmds.append(new_cmd1 + name + new_cmd2 +
                    name + new_cmd3 +
                    label + new_cmd4)
            return cmds

        def make_edges(edges):
            all_paths = []
            for edge in edges:
                path = "(" + edge[0] + ")"
                for i in range(1, len(edge)):
                    path += " edge (%s)" % edge[i]
                all_paths.append(path)
            return all_paths
        return ("{ " +
            "".join(make_cmd(nodes)) +
            begin_env +
                (matrix_begin +
                    "\\\\ \n".join(matrix) +
                matrix_end +
                ("\n" +
                path_begin +
                    "\n\t".join(make_edges(edges)) +
                path_end if len(edges) > 0 else "")
                if len(matrix) > 0 else "") +
            end_env +
            "}")
