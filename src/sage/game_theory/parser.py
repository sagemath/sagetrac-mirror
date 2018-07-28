from __future__ import print_function

class Parser():
    r"""
    A class for parsing the outputs of different algorithms called in other
    software packages.

    Two parsers are included, one for the ``'lrs'`` algorithm and another for
    the ``'LCP'`` algorithm.
    """

    def __init__(self, raw_string):
        """
        Initialise a Parser instance by storing a raw_string
        (currently only used with H representation of a game).

        TESTS:

        Simply checking that we have the correct string output
        for the H representation (which is the format required
        for the ``'lrs'`` algorithm)::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: raw_string = g._Hrepresentation(A, -A)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 1 0 0
            0 0 1 0
            0 1 3  1
            0 2 2  1
            -1 1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 -1 -2  1
            0 -3 -2  1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
            <BLANKLINE>

        The specific case of a game with 1 strategy for each player::

            sage: A = matrix([[1]])
            sage: B = matrix([[5]])
            sage: g = NormalFormGame([A,B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 1 0
            0 -5  1
            -1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 -1  1
            0 1 0
            -1 1 0
            end
            <BLANKLINE>

        Another test::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0  9 -6  4  1
            0 -7  2 -6  1
            0 -9  3  10  1
            -1 1 1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 7 5 -5  1
            0 -5 -5 -3  1
            0 -1 6 -1  1
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            -1 1 1 1 0
            end
            <BLANKLINE>

        This class is also used to parse the output of algorithms from the gambit
        python interface using the `format_gambit` function.
        """
        self.raw_string = raw_string

    def format_lrs(self):
        """
        Parses the output of lrs so as to return vectors
        corresponding to equilibria.

        TESTS::

            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: game1_str, game2_str = g._Hrepresentation(A, -A)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = open(g1_name, 'w')
            sage: g2_file = open(g2_name, 'w')
            sage: _ = g1_file.write(game1_str)
            sage: g1_file.close()
            sage: _ = g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['lrsnash', g1_name, g2_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [row for row in process.stdout]  # optional - lrslib

        The above creates a game, writes the H representation to
        temporary files, calls lrs and stores the output in `lrs_output`
        (here slicing to get rid of some system parameters that get returned)::

            sage: lrs_output[5:16]  # optional - lrslib
            ['\n',
             '***** 4 4 rational\n',
             '2  0  1  2 \n',
             '1  1/2  1/2 -2 \n',
             '\n',
             '2  0  1  2 \n',
             '1  0  1 -2 \n',
             '\n',
             '\n',
             '*Number of equilibria found: 2\n',
             '*Player 1: vertices=3 bases=3 pivots=5\n']

        The above is pretty messy, here is the output when we put it through
        the parser::

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrslib
            sage: nasheq  # optional - lrslib
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]

        Another game::

            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: game1_str, game2_str = g._Hrepresentation(A, B)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = open(g1_name, 'w')
            sage: g2_file = open(g2_name, 'w')
            sage: _ = g1_file.write(game1_str)
            sage: g1_file.close()
            sage: _ = g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['lrsnash', g1_name, g2_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [row for row in process.stdout]  # optional - lrslib
            sage: print(lrs_output[5:20])  # optional - lrslib
            ['\n',
             '***** 5 5 rational\n',
             '2  1/7  0  6/7  23/7 \n',
             '2  0  1/6  5/6  10/3 \n',
             '1  1/3  2/3  0  1 \n',
             '\n',
             '2  0  0  1  5 \n',
             '1  1  0  0  9 \n',
             '\n',
             '2  1  0  0  5 \n',
             '1  0  1  0  6 \n',
             '\n',
             '\n',
             '*Number of equilibria found: 4\n',
             '*Player 1: vertices=6 bases=7 pivots=10\n']

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrslib
            sage: sorted(nasheq)  # optional - lrslib
            [[(0, 1, 0), (1, 0, 0)],
             [(1/3, 2/3, 0), (0, 1/6, 5/6)],
             [(1/3, 2/3, 0), (1/7, 0, 6/7)],
             [(1, 0, 0), (0, 0, 1)]]
        """
        equilibria = []
        from sage.misc.sage_eval import sage_eval
        from itertools import groupby
        for collection in [list(x[1]) for x in groupby(self.raw_string[7:], lambda x: x == '\n')]:
            if collection[0].startswith('2'):
                s1 = tuple([sage_eval(k) for k in collection[-1].split()][1:-1])
                for s2 in collection[:-1]:
                    s2 = tuple([sage_eval(k) for k in s2.split()][1:-1])
                    equilibria.append([s1, s2])

        return equilibria

    def format_gambit(self, gambit_game):
        """
        Parses the output of gambit so as to return vectors
        corresponding to equilibria obtained using the LCP algorithm.

        TESTS:

        Here we construct a two by two game in gambit::

            sage: import gambit  # optional - gambit
            sage: from sage.game_theory.parser import Parser
            sage: g = gambit.Game.new_table([2,2])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(2)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(1)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(0)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(0)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(0)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(0)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(2)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the output of the LCP algorithm::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0], [1.0, 0.0]]>,
             <NashProfile for '': [[0.6666666667, 0.3333333333], [0.3333333333, 0.6666666667]]>,
             <NashProfile for '': [[0.0, 1.0], [0.0, 1.0]]>]

        The Parser class outputs the equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)], [(0.6666666667, 0.3333333333), (0.3333333333, 0.6666666667)], [(0.0, 1.0), (0.0, 1.0)]]

        Here is another game::

            sage: g = gambit.Game.new_table([2,2])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(4)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(8)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(0)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(1)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(3)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(0)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the LCP output::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0], [1.0, 0.0]]>]

        The corresponding parsed equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)]]

        Here is a larger degenerate game::

            sage: g = gambit.Game.new_table([3,3])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(-7)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(-9)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(-5)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(7)  # optional - gambit
            sage: g[int(0), int(2)][int(0)] = int(5)  # optional - gambit
            sage: g[int(0), int(2)][int(1)] = int(9)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(5)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(6)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(5)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(-2)  # optional - gambit
            sage: g[int(1), int(2)][int(0)] = int(3)  # optional - gambit
            sage: g[int(1), int(2)][int(1)] = int(-3)  # optional - gambit
            sage: g[int(2), int(0)][int(0)] = int(1)  # optional - gambit
            sage: g[int(2), int(0)][int(1)] = int(-4)  # optional - gambit
            sage: g[int(2), int(1)][int(0)] = int(-6)  # optional - gambit
            sage: g[int(2), int(1)][int(1)] = int(6)  # optional - gambit
            sage: g[int(2), int(2)][int(0)] = int(1)  # optional - gambit
            sage: g[int(2), int(2)][int(1)] = int(-10)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the LCP output::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]>,
             <NashProfile for '': [[0.3333333333, 0.6666666667, 0.0], [0.1428571429, 0.0, 0.8571428571]]>,
             <NashProfile for '': [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]>]

        The corresponding parsed equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
             [(0.3333333333, 0.6666666667, 0.0), (0.1428571429, 0.0, 0.8571428571)],
             [(0.0, 1.0, 0.0), (1.0, 0.0, 0.0)]]

        Note, that this differs from the same output of the lrs algorithm due
        the fact that the game is degenerate.
        """
        nice_stuff = []
        for gambitstrategy in self.raw_string:  #looks at each individual profile
            gambitstrategy = list(gambitstrategy)
            profile = [tuple(gambitstrategy[:len(gambit_game.players[int(0)].strategies)])]
            for player in list(gambit_game.players)[1:]:
                previousplayerstrategylength = len(profile[-1])
                profile.append(tuple(gambitstrategy[previousplayerstrategylength: previousplayerstrategylength + len(player.strategies)]))
            nice_stuff.append(profile)

        return nice_stuff

    def format_gambit_efg_tree(self, gambit_game):
        """
        Parses the output of gambit so as to returns a list of each equilibria,
        with tuples that contains information sets and corresponding dictionaries
        that maps actions to probabilities, corresponding to equilibria
        obtained using Gambit's LCP algorithm.

        TESTS:

        Here we construct a two person extensive form game in Gambit::

            sage: import gambit  # optional - gambit
            sage: from sage.game_theory.parser import Parser
            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: g.title = "Parser example"  # optional - gambit
            sage: iset = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: iset.label = "a"  # optional - gambit
            sage: iset.actions[int(0)].label = "X"  # optional - gambit
            sage: iset.actions[int(1)].label = "W"  # optional - gambit
            sage: iset = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: iset.label = "b"  # optional - gambit
            sage: iset.actions[int(0)].label = "D"  # optional - gambit
            sage: iset.actions[int(1)].label = "C"  # optional - gambit
            sage: iset = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: iset.label = "c"  # optional - gambit
            sage: iset.actions[int(0)].label = "B"  # optional - gambit
            sage: iset.actions[int(1)].label = "A"  # optional - gambit
            sage: iset = g.root.children[int(0)].children[int(1)].append_move(g.players["1"], int(2))  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].label = 'Node 3'  # optional - gambit
            sage: iset.label = "d"  # optional - gambit
            sage: iset.actions[int(0)].label = "Z"  # optional - gambit
            sage: iset.actions[int(1)].label = "Y"  # optional - gambit
            sage: XD = g.outcomes.add("XD")  # optional - gambit
            sage: XD[int(0)] = int(2)  # optional - gambit
            sage: XD[int(1)] = int(0)  # optional - gambit
            sage: XCZ = g.outcomes.add("XCZ")  # optional - gambit
            sage: XCZ[int(0)] = int(3)  # optional - gambit
            sage: XCZ[int(1)] = int(1)  # optional - gambit
            sage: XCY = g.outcomes.add("XCY")  # optional - gambit
            sage: XCY[int(0)] = int(4)  # optional - gambit
            sage: XCY[int(1)] = int(2)  # optional - gambit
            sage: WB = g.outcomes.add("WB")  # optional - gambit
            sage: WB[int(0)] = int(3)  # optional - gambit
            sage: WB[int(1)] = int(5)  # optional - gambit
            sage: WA = g.outcomes.add("WA")  # optional - gambit
            sage: WA[int(0)] = int(4)  # optional - gambit
            sage: WA[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(0)].outcome = XCZ  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(1)].outcome = XCY  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = WB  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = WA  # optional - gambit

        Here is the output of the LCP algorithm::

            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit  # optional - gambit
            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit  # optional - gambit
            [<NashProfile for 'Parser example': [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]>]

        The Parser class outputs the equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit_efg_tree(g)  # optional - gambit
            sage: expected_outcome = [((('Root',), {'W': 0.0, 'X': 1.0}), (('Node 3',), {'Y': 1.0, 'Z': 0.0}),
            ....: (('Node 1',), {'C': 0.0, 'D': 1.0}), (('Node 2',), {'A': 1.0, 'B': 0.0}))]
            sage: nasheq == expected_outcome # optional - gambit
            True

        If we change one of the outputs for the above tree, more nash equilibria are obtained::

            sage: alternate_output = g.outcomes.add()  # optional - gambit
            sage: alternate_output[int(0)] = int(5)  # optional - gambit
            sage: alternate_output[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = alternate_output  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit
            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for 'Parser example': [0.0, 1.0, 0.5, 0.5, 0.0, 1.0, 1.0, 0.0]>,
             <NashProfile for 'Parser example': [0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0]>,
             <NashProfile for 'Parser example': [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]>]
            sage: nasheq = Parser(LCP_output).format_gambit_efg_tree(g)  # optional - gambit
            sage: expected_outcome = [((('Root',), {'W': 1.0, 'X': 0.0}), (('Node 3',), {'Y': 0.5, 'Z': 0.5}),
            ....: (('Node 1',), {'C': 1.0, 'D': 0.0}), (('Node 2',), {'A': 0.0, 'B': 1.0})),
            ....: ((('Root',), {'W': 1.0, 'X': 0.0}), (('Node 3',), {'Y': 0.5, 'Z': 0.5}),
            ....: (('Node 1',), {'C': 0.5, 'D': 0.5}), (('Node 2',), {'A': 0.0, 'B': 1.0})),
            ....: ((('Root',), {'W': 0.0, 'X': 1.0}), (('Node 3',), {'Y': 1.0, 'Z': 0.0}),
            ....: (('Node 1',), {'C': 0.0, 'D': 1.0}), (('Node 2',), {'A': 1.0, 'B': 0.0}))]
            sage: nasheq  == expected_outcome  # optional - gambit
            True

        Another test::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: iset = g.root.append_move(g.players["1"], int(3))  # optional - gambit
            sage: g.root.label = 'Root' # optional - gambit
            sage: iset.actions[int(0)].label = "A"  # optional - gambit
            sage: iset.actions[int(1)].label = "B"  # optional - gambit
            sage: iset.actions[int(2)].label = "C"  # optional - gambit
            sage: iset = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: iset.actions[int(0)].label = "D"  # optional - gambit
            sage: iset.actions[int(1)].label = "E"  # optional - gambit
            sage: iset = g.root.children[int(0)].children[int(1)].append_move(g.players["1"], int(2))  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].label = 'Node 3' # optional - gambit
            sage: iset.actions[int(0)].label = "F"  # optional - gambit
            sage: iset.actions[int(1)].label = "G"  # optional - gambit
            sage: iset = g.root.children[int(2)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(2)].label = 'Node 2'  # optional - gambit
            sage: iset.actions[int(0)].label = "H"  # optional - gambit
            sage: iset.actions[int(1)].label = "I"  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(9)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(2)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(2)].children[int(1)].outcome = outcome  # optional - gambit

        The output of the LCP algorithm::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [0.0, 1.0, 0.0, 0.5, 0.5, 0.75, 0.25, 0.0, 1.0]>]

        The output of the Parser::
            sage: nasheq = Parser(LCP_output).format_gambit_efg_tree(g)  # optional - gambit
            sage: expected_outcome = [((('Root',), {'A': 0.0, 'B': 1.0, 'C': 0.0}),
            ....: (('Node 3',), {'F': 0.5, 'G': 0.5}),
            ....: (('Node 1',), {'D': 0.75, 'E': 0.25}), (('Node 2',), {'H': 0.0, 'I': 1.0}))]
            sage: nasheq == expected_outcome  # optional - gambit
            True

        Another test with a different tree::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: iset = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: iset.actions[int(0)].label = "A"  # optional - gambit
            sage: iset.actions[int(1)].label = "B"  # optional - gambit
            sage: iset = g.root.children[int(0)].append_move(g.players["2"], int(3))  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: iset.actions[int(0)].label = "C"  # optional - gambit
            sage: iset.actions[int(1)].label = "D"  # optional - gambit
            sage: iset.actions[int(2)].label = "E"  # optional - gambit
            sage: iset = g.root.children[int(1)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: iset.actions[int(0)].label = "F"  # optional - gambit
            sage: iset.actions[int(1)].label = "G"  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(9)  # optional - gambit
            sage: outcome[int(1)] = int(1)  # optional - gambit
            sage: g.root.children[int(0)].children[int(2)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        The output of the LCP algorithm::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [0.0, 1.0, 0.75, 0.25, 0.0, 0.0, 1.0]>]

        The output of the Parser::

            sage: nasheq = Parser(LCP_output).format_gambit_efg_tree(g)  # optional - gambit
            sage: expected_outcome = [((('Root',), {'A': 0.0, 'B': 1.0}),
            ....: (('Node 1',), {'C': 0.75, 'D': 0.25, 'E': 0.0}),
            ....: (('Node 2',), {'F': 0.0, 'G': 1.0}))]
            sage: nasheq == expected_outcome  # optional - gambit
            True

        Another test::

            sage: g = gambit.Game.new_tree()  # optional - gambit
            sage: g.players.add("1")  # optional - gambit
            <Player [0] '1' in game ''>
            sage: g.players.add("2")  # optional - gambit
            <Player [1] '2' in game ''>
            sage: iset = g.root.append_move(g.players["1"], int(2))  # optional - gambit
            sage: g.root.label = 'Root'  # optional - gambit
            sage: iset.actions[int(0)].label = "A"  # optional - gambit
            sage: iset.actions[int(1)].label = "B"  # optional - gambit
            sage: iset = g.root.children[int(0)].append_move(g.players["2"], int(2))  # optional - gambit
            sage: g.root.children[int(0)].label = 'Node 1'  # optional - gambit
            sage: iset.actions[int(0)].label = "C"  # optional - gambit
            sage: iset.actions[int(1)].label = "D"  # optional - gambit
            sage: g.root.children[int(1)].append_move(iset)  # optional - gambit
            <Infoset [0] '' for player '2' in game ''>
            sage: g.root.children[int(1)].label = 'Node 2'  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(1)  # optional - gambit
            sage: outcome[int(1)] = int(5)  # optional - gambit
            sage: g.root.children[int(0)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(5)  # optional - gambit
            sage: outcome[int(1)] = int(2)  # optional - gambit
            sage: g.root.children[int(0)].children[int(1)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(3)  # optional - gambit
            sage: outcome[int(1)] = int(0)  # optional - gambit
            sage: g.root.children[int(1)].children[int(0)].outcome = outcome  # optional - gambit
            sage: outcome = g.outcomes.add()  # optional - gambit
            sage: outcome[int(0)] = int(2)  # optional - gambit
            sage: outcome[int(1)] = int(7)  # optional - gambit
            sage: g.root.children[int(1)].children[int(1)].outcome = outcome  # optional - gambit

        The output of the LCP algorithm::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [0.7, 0.3, 0.6, 0.4]>]

        The output of the Parser::

            sage: nasheq = Parser(LCP_output).format_gambit_efg_tree(g)  # optional - gambit
            sage: expected_outcome = [((('Root',), {'A': 0.7, 'B': 0.3}),
            ....: (('Node 1', 'Node 2',), {'C': 0.6, 'D': 0.4}))]
            sage: nasheq  == expected_outcome # optional - gambit
            True
        """

        nice_stuff = []
        for gambitstrategy in self.raw_string:
            gambitstrategy = list(gambitstrategy)
            infoset_action_count = int(0)
            gambitstrategy_list = []
            for player in list(gambit_game.players):
                for infoset in list(player.infosets):
                    infoset_strategy = gambitstrategy[infoset_action_count: infoset_action_count + int(len(list(infoset.actions)))]
                    infoset_action_count += int(len(infoset.actions))
                    action_dict = {action.label:infoset_strategy[i] for i,action in enumerate(infoset.actions)}
                    node_list = []
                    for node in infoset.members:
                        node_list.append(node.label)
                    gambitstrategy_list.append((tuple(node_list), action_dict))
            nice_stuff.append(tuple(gambitstrategy_list))
        return nice_stuff
