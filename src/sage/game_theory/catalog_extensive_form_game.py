r"""
A catalog of extensive form games.

This allows us to construct common games directly::

    sage: g = game_theory.extensive_form_game.PrisonersDilemma()
    sage: g
    Prisoners Dilemma - Extensive Form Game with the following underlying tree: {...}

We can then immediately obtain the Nash equilibrium for this game::

    sage: expected_nash = [((('a',), {'Defect': 1.0, 'Cooperate': 0.0}), (('b', 'c'), {'Defect': 1.0, 'Cooperate': 0.0}))]
    sage: expected_nash == g.obtain_nash() # optional - gambit
    True

We can also look at all the outputs and attributes of the game::

    sage: g.nodes
    [EFG Node "a", EFG Node "b", EFG Node "c"]
    sage: g.leafs
    [EFG Leaf with utilities - (-4, -4),
     EFG Leaf with utilities - (0, -5),
     EFG Leaf with utilities - (-5, 0),
     EFG Leaf with utilities - (-2, -2)]
    sage: g.tree_root
    EFG Node "a"
    sage: g.players
    [EFG Player "Prisoner one", EFG Player "Prisoner two"]
    sage: g.info_sets
    [[EFG Node "a"], [EFG Node "b", EFG Node "c"]]

The docstrings give an interpretation of each game.

More information is available in the following references:

REFERENCES:

.. [Basu] Kaushik Basu.
   *The Traveler's Dilemma: Paradoxes of Rationality in Game Theory*.
   The American Economic Review (1994): 391-395.

.. [Cressman] Cressman, Ross.
   *Evolutionary dynamics and extensive form games*.
   MIT Press.

.. [McMillan] John McMillan.
   *Games, strategies, and managers*.
   Oxford University Press.

.. [Skyrms] Brian Skyrms.
   *The stag hunt and the evolution of social structure*.
   Cambridge University Press.

.. [Watson] Joel Watson.
   *Strategy: an introduction to game theory*.
   WW Norton.

.. [Webb] James Webb.
   *Game theory: decisions, interaction and Evolution*.
   Springer Science & Business Media.

AUTHOR:

- Francis Rhys Ward and Vince Knight (06-2014)
"""
from sage.game_theory.extensive_form_game import ExtensiveFormGame


def PrisonersDilemma(R=-2, P=-4, S=-5, T=0):
    r"""
    Return a Prisoners dilemma game.

    Assume two thieves have been caught by the police
    and separated for questioning.
    If both thieves cooperate and do not divulge any information they will
    each get a short sentence.
    If one defects he/she is offered a deal while the other thief will get a
    long sentence.
    If they both defect they both get a medium length sentence. [Webb]_

    This can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.PrisonersDilemma()
        p = Prisoners_Dilema.plot(view_info_sets=True)
        sphinx_plot(p)

    Where `T > R > P > S`.

    - `R` denotes the reward received for cooperating.
    - `S` denotes the 'sucker' utility.
    - `P` denotes the utility for punishing the other player.
    - `T` denotes the temptation payoff.
    
    This is the default game created::

        sage: g = game_theory.extensive_form_game.PrisonersDilemma()
        sage: g
        Prisoners Dilemma - Extensive Form Game with the following underlying tree: {...}

    When we test whether the functions are returning the games that we expect,
    we will be testing a dictionary whith keys as nodes, and values as the nodes' children names.
    In this way we will test that the tree created by the function is correct::

        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        
    We can also test the utilities at each leaf::
    
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[-4, -4], [0, -5], [-5, 0], [-2, -2]]

    There is a single Nash equilibrium for this at which both thieves defect.
    This can be implemented in Sage using the following::

        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'Defect': 1.0, 'Cooperate': 0.0}), (('b', 'c'), {'Defect': 1.0, 'Cooperate': 0.0}))]
        sage: N == d
        True

    Note that we can pass other values of R, P, S, T::

        sage: g = game_theory.extensive_form_game.PrisonersDilemma(R=-1, P=-2, S=-3, T=0)
        sage: g
        Prisoners Dilemma - Extensive Form Game with the following underlying tree: {...}
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[-2, -2], [0, -3], [-3, 0], [-1, -1]]

    If we pass values that fail the defining requirement: `T > R > P > S`
    we get an error message::

        sage: g = game_theory.extensive_form_game.PrisonersDilemma(R=-1, P=-2, S=0, T=5)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a Prisoners Dilemma must be
        of the form T > R > P > S

    """
    if not (T > R > P > S):
        raise TypeError("the input values for a Prisoners Dilemma must be of the form T > R > P > S")
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Prisoner one')
    player_2 = EFG_Player('Prisoner two')
    leaf_1 = EFG_Leaf({player_1: R, player_2: R})
    leaf_2 = EFG_Leaf({player_1: S, player_2: T})
    leaf_3 = EFG_Leaf({player_1: T, player_2: S})
    leaf_4 = EFG_Leaf({player_1: P, player_2: P})
    node_3 = EFG_Node(player_2, {'Cooperate': leaf_1, 'Defect': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'Cooperate': leaf_3, 'Defect': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'Cooperate': node_3, 'Defect': node_2}, 'a')
    g = Prisoners_Dilema = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3])
    g.rename('Prisoners Dilemma - ' + repr(g))
    return g


def CoordinationGame(A=10, a=5, B=0, b=0, C=0, c=0, D=5, d=10):
    r"""
    Return a Coordination Game with two players whom each have 2 strategies.

    A coordination game is a particular type of game where the pure Nash
    equilibrium is for the players to pick the same strategies [Webb]_.

    An often used example of this can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.CoordinationGame()
        p = CoordinationGame.plot(view_info_sets=True)
        sphinx_plot(p)

    Where `A > B, D > C` and `a > c, d > b`.

    This is the default version of the game created by this function::

        sage: g = game_theory.extensive_form_game.CoordinationGame()
        sage: g
        Coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[5, 10], [0, 0], [0, 0], [10, 5]]
        
    There are two pure Nash equilibria and one mixed::

        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'Y': 0.0, 'X': 1.0}), (('b', 'c'), {'Y': 0.0, 'X': 1.0})), ((('a',), {'Y': 0.3333333333, 'X': 0.6666666667}), (('b', 'c'), {'Y': 0.6666666667, 'X': 0.3333333333})), ((('a',), {'Y': 1.0, 'X': 0.0}), (('b', 'c'), {'Y': 1.0, 'X': 0.0}))]
        sage: N == d
        True

    We can also pass different values of the input parameters::

        sage: g = game_theory.extensive_form_game.CoordinationGame(A=9, a=6,
        ....:                                B=2, b=1, C=0, c=1, D=4, d=11)
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[4, 11], [2, 1], [0, 1], [9, 6]]

    Note that an error is returned if the defining inequalities are
    not obeyed `A > B, D > C` and `a > c, d > b`::

        sage: g = game_theory.extensive_form_game.CoordinationGame(A=9, a=6,
        ....:                               B=0, b=1, C=2, c=10, D=4, d=11)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a Coordination game must
                        be of the form A > B, D > C, a > c and d > b
    """
    if not (A > B  and  D > C and a > c and d > b):
        raise TypeError("the input values for a Coordination game must be of the form A > B, D > C, a > c and d > b")
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf
    
    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: A, player_2: a})
    leaf_2 = EFG_Leaf({player_1: B, player_2: b})
    leaf_3 = EFG_Leaf({player_1: C, player_2: c})
    leaf_4 = EFG_Leaf({player_1: D, player_2: d})
    node_3 = EFG_Node(player_2, {'X': leaf_1, 'Y': leaf_3}, 'c')
    node_2 = EFG_Node(player_2, {'X': leaf_2, 'Y': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'X': node_3, 'Y': node_2}, 'a')
    g = Coordination_game = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3])
    g.rename('Coordination game - ' + repr(g))
    return g


def BattleOfTheSexes():
    r"""
    Return a Battle of the Sexes game.

    Consider two payers: Amy and Bob.
    Amy prefers to play video games and Bob prefers to
    watch a movie. They both however want to spend their evening
    together. [Webb]_

    This can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.BattleOfTheSexes()
        p = BattleOfTheSexes.plot(view_info_sets=True)
        sphinx_plot(p)

    This is a particular type of Coordination Game.
    There are three Nash equilibria:

    1. Amy and Bob both play video games;
    2. Amy and Bob both watch a movie;
    3. Amy plays video games 75% of the time and Bob watches a movie 75% of the time.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.BattleOfTheSexes()
        sage: g
        Battle of the sexes - Coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'Y': 0.0, 'X': 1.0}), (('b', 'c'), {'Y': 0.0, 'X': 1.0})), ((('a',), {'Y': 0.25, 'X': 0.75}), (('b', 'c'), {'Y': 0.75, 'X': 0.25})), ((('a',), {'Y': 1.0, 'X': 0.0}), (('b', 'c'), {'Y': 1.0, 'X': 0.0}))]
        sage: N == d
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[2, 3], [0, 0], [1, 1], [3, 2]]
    """
    g = CoordinationGame(A=3, a=2, B=0, b=0, C=1, c=1, D=2, d=3)
    g.rename('Battle of the sexes - ' + repr(g))
    return g


def StagHunt():
    r"""
    Return a Stag Hunt game.

    Assume two friends go out on a hunt. Each can individually choose to hunt
    a stag or hunt a hare. Each player must choose an action without knowing
    the choice of the other. If an individual hunts a stag, he must have the
    cooperation of his partner in order to succeed. An individual can get a
    hare by himself, but a hare is worth less than a stag. [Skyrms]_

    This can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.StagHunt()
        p = StagHunt.plot(view_info_sets=True)
        sphinx_plot(p)

    This is a particular type of Coordination Game.
    There are three Nash equilibria:

        1. Both friends hunting the stag.
        2. Both friends hunting the hare.
        3. Both friends hunting the stag 2/3rds of the time.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.StagHunt()
        sage: g
        Stag hunt - Coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'X': 1.0, 'Y': 0.0}), (('b', 'c'), {'X': 1.0, 'Y': 0.0})), ((('a',), {'X': 0.6666666667, 'Y': 0.3333333333}), (('b', 'c'), {'X': 0.6666666667, 'Y': 0.3333333333})), ((('a',), {'X': 0.0, 'Y': 1.0}), (('b', 'c'), {'X': 0.0, 'Y': 1.0}))]
        sage: N == d
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[2, 2], [4, 0], [0, 4], [5, 5]]


    """
    g = CoordinationGame(A=5, a=5, B=4, b=0, C=0, c=4, D=2, d=2)
    g.rename('Stag hunt - ' + repr(g))
    return g


def AntiCoordinationGame(A=3, a=3, B=5, b=1, C=1, c=5, D=0, d=0):
    r"""
    Return an AntiCoordination Game with two players whom each have 2 strategies.

    An anti coordination game is a particular type of game where the pure Nash
    equilibria is for the players to pick different strategies.

    The default example of this can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.AntiCoordinationGame()
        p = AntiCoordinationGame.plot(view_info_sets=True)
        sphinx_plot(p)

    Where `A < B, D < C` and `a < c, d < b`.

    This is the default version of the game created::

        sage: g = game_theory.extensive_form_game.AntiCoordinationGame()
        sage: g
        Anti coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[0, 0], [5, 1], [1, 5], [3, 3]]

    There are two pure Nash equilibria and one mixed::

        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'Y': 0.0, 'X': 1.0}), (('b', 'c'), {'Y': 1.0, 'X': 0.0})), ((('a',), {'Y': 0.6666666667, 'X': 0.3333333333}), (('b', 'c'), {'Y': 0.6666666667, 'X': 0.3333333333})), ((('a',), {'Y': 1.0, 'X': 0.0}), (('b', 'c'), {'Y': 0.0, 'X': 1.0}))]
        sage: N == d
        True

    We can also pass different values of the input parameters::

        sage: g = game_theory.extensive_form_game.AntiCoordinationGame(A=1, a=1,
        ....:                                     B=2, b=8, C=4, c=2, D=1, d=0)
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[1, 0], [2, 8], [4, 2], [1, 1]]

    Note that an error is returned if the defining inequality is
    not obeyed `A > B, D > C` and `a > c, d > b`::

        sage: g = game_theory.extensive_form_game.AntiCoordinationGame(A=8, a=3,
        ....:                                     B=4, b=2, C=2, c=8, D=1, d=0)
        Traceback (most recent call last):
        ...
        TypeError: the input values for an Anti coordination game must be of the form A < B, D < C, a < c and d < b
    """
    if not (A < B  and  D < C and a < c and d < b):
        raise TypeError("the input values for an Anti coordination game must be of the form A < B, D < C, a < c and d < b")
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: A, player_2: a})
    leaf_2 = EFG_Leaf({player_1: B, player_2: b})
    leaf_3 = EFG_Leaf({player_1: C, player_2: c})
    leaf_4 = EFG_Leaf({player_1: D, player_2: d})
    node_3 = EFG_Node(player_2, {'X': leaf_1, 'Y': leaf_3}, 'c')
    node_2 = EFG_Node(player_2, {'X': leaf_2, 'Y': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'X': node_3, 'Y': node_2}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3])
    g.rename('Anti coordination game - ' + repr(g))
    return g


def HawkDove(v=2, c=3):
    r"""
    Return a Hawk Dove game.

    Suppose two birds of prey must share a limited resource `v`.
    The birds can act like a hawk or a dove.

    - If a dove meets a hawk, the hawk takes the resources.
    - If two doves meet they share the resources.
    - If two hawks meet, one will win (with equal expectation) and take the
      resources while the other will suffer a cost of `c` where
      `c>v`. [Webb]_

    The general form of this game can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.HawkDove(A=v/2-c, a=v/2-c, B=0, b=v, C=v, c=0, D=v/2, d=v/2)
        p = HawkDove.plot(view_info_sets=True)
        sphinx_plot(p)

    Here is the default game with the values of `v=2` and `c=3`.

        sage: g = game_theory.extensive_form_game.HawkDove()
        sage: g
        Hawk-Dove - Anti coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[1, 1], [0, 2], [2, 0], [-2, -2]]

    This is a particular example of an anti coordination game.
    There are three Nash equilibria:

        1. One bird acts like a Hawk and the other like a Dove.
        2. Both birds mix being a Hawk and a Dove

    This can be implemented in Sage using the following::

        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'X': 1.0, 'Y': 0.0}), (('b', 'c'), {'X': 0.0, 'Y': 1.0})), ((('a',), {'X': 0.3333333333, 'Y': 0.6666666667}), (('b', 'c'), {'X': 0.3333333333, 'Y': 0.6666666667})), ((('a',), {'X': 0.0, 'Y': 1.0}), (('b', 'c'), {'X': 1.0, 'Y': 0.0}))]
        sage: N == d
        True

    Note that an error is returned if the defining inequality is not obeyed
    `c < v`:

        sage: g = game_theory.extensive_form_game.HawkDove(v=5, c=1)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a Hawk Dove game must be of the form c > v
    """
    if not (c>v):
        raise TypeError("the input values for a Hawk Dove game must be of the form c > v")
    g = AntiCoordinationGame(A=v/2-c, a=v/2-c, B=0, b=v, C=v, c=0, D=v/2, d=v/2)
    g.rename('Hawk-Dove - ' + repr(g))
    return g


def Pigs():
    r"""
    Return a Pigs game.

    Consider two pigs.
    One dominant pig and one subservient pig.
    These pigs share a pen.
    There is a lever in the pen that delivers 6 units of food but if either pig
    pushes the lever it will take them a little while to get to the food as well
    as cost them 1 unit of food.
    If the dominant pig pushes the lever, the subservient pig has some time
    to eat two thirds of the food before being pushed out of the way.
    If the subservient pig pushes the lever,
    the dominant pig will eat all the food.
    Finally if both pigs go to push the lever the subservient pig will be able
    to eat a third of the food (and they will also both lose 1 unit of food).
    [McMillan]_

    The default form of this game can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.Pigs()
        p = Pigs.plot(view_info_sets=True)
        sphinx_plot(p)

    There is a single Nash equilibrium at which the dominant pig pushes the
    lever and the subservient pig does not.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.Pigs()
        sage: g
        Pigs - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[0, 0], [-1, 6], [4, 1], [1, 3]]
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'DPL': 0.0, 'PL': 1.0}), (('b', 'c'), {'DPL': 1.0, 'PL': 0.0}))]
        sage: N == d
        True
    """
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf
    
    player_1 = EFG_Player('Dominant pig')
    player_2 = EFG_Player('Subserviant pig')
    leaf_1 = EFG_Leaf({player_1: 3, player_2: 1})
    leaf_2 = EFG_Leaf({player_1: 1, player_2: 4})
    leaf_3 = EFG_Leaf({player_1: 6, player_2: -1})
    leaf_4 = EFG_Leaf({player_1: 0, player_2: 0})
    node_3 = EFG_Node(player_2, {'PL': leaf_1, 'DPL': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'PL': leaf_3, 'DPL': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'PL': node_3, 'DPL': node_2}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3])
    g.rename('Pigs - ' + repr(g))
    return g


def MatchingPennies():
    r"""
    Return a Matching Pennies game.

    Consider two players who can choose to display a coin either Heads
    facing up or Tails facing up.
    If both players show the same face then player 1 wins,
    if not then player 2 wins. [Webb]_

    This can be modeled as a zero extensive form game with the following
    tree:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.MatchingPennies()
        p = MatchingPennies.plot(view_info_sets=True)
        sphinx_plot(p)

    There is a single Nash equilibria at which both players randomly
    (with equal probability) pick heads or tails.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.MatchingPennies()
        sage: g
        Matching pennies - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['c', 'b'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[-1, 1], [1, -1], [1, -1], [-1, 1]]
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'H': 0.5, 'T': 0.5}), (('b', 'c'), {'H': 0.5, 'T': 0.5}))]
        sage: N == d
        True
    """
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: 1, player_2: -1})
    leaf_2 = EFG_Leaf({player_1: -1, player_2: 1})
    leaf_3 = EFG_Leaf({player_1: -1, player_2: 1})
    leaf_4 = EFG_Leaf({player_1: 1, player_2: -1})
    node_3 = EFG_Node(player_2, {'H': leaf_1, 'T': leaf_2}, 'c')
    node_2 = EFG_Node(player_2, {'H': leaf_3, 'T': leaf_4}, 'b')
    node_1 = EFG_Node(player_1, {'H': node_3, 'T': node_2}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3])
    g.rename('Matching pennies - ' + repr(g))
    return g


def RPS():
    r"""
    Return a Rock-Paper-Scissors game.

    Rock-Paper-Scissors is a zero sum game usually played between two
    players where each player simultaneously forms one of three
    shapes with an outstretched hand.The game has only three possible outcomes
    other than a tie: a player who decides to play rock will beat another
    player who has chosen scissors ("rock crushes scissors") but will lose to
    one who has played paper ("paper covers rock"); a play of paper will lose
    to a play of scissors ("scissors cut paper"). If both players throw the
    same shape, the game is tied and is usually immediately replayed to break
    the tie. [Webb]_

    This can be modeled as a zero sum extensive form game with the following
    tree:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.RPS()
        p = RPS.plot(view_info_sets=True)
        sphinx_plot(p)

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.RPS()
        sage: g
        Rock-Paper-Scissors - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['c', 'd', 'b'], 'b': ['Leaf 1', 'Leaf 2', 'Leaf 3'], 'c': ['Leaf 4', 'Leaf 5', 'Leaf 6'], 'd': ['Leaf 7', 'Leaf 8', 'Leaf 9']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[0, 1], [1, 0], [0, 0], [0, 0], [0, 1], [1, 0], [1, 0], [0, 0], [0, 1]]
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'P': 0.3333333333, 'S': 0.3333333333, 'R': 0.3333333333}), (('c', 'b', 'd'), {'P': 0.3333333333, 'S': 0.3333333333, 'R': 0.3333333333}))]
        sage: N == d
        True
    """
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_4 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_5 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_6 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_7 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_8 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_9 = EFG_Leaf({player_1: 0, player_2: 0})
    node_4 = EFG_Node(player_2, {'R': leaf_7, 'P': leaf_8, 'S': leaf_9 }, 'd')
    node_3 = EFG_Node(player_2, {'R': leaf_4, 'P': leaf_5, 'S': leaf_6 }, 'c')
    node_2 = EFG_Node(player_2, {'R': leaf_1, 'P': leaf_2, 'S': leaf_3}, 'b')
    node_1 = EFG_Node(player_1, {'R': node_2, 'P': node_3, 'S': node_4}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3, node_4])
    g.rename('Rock-Paper-Scissors - ' + repr(g))
    return g

def RPSLSp():
    r"""
    Return a Rock-Paper-Scissors-Lizard-Spock game.

    `Rock-Paper-Scissors-Lizard-Spock
    <http://www.samkass.com/theories/RPSSL.html>`_ is an extension of
    Rock-Paper-Scissors.
    It is a zero sum game usually played between two
    players where each player simultaneously forms one of five
    shapes with an outstretched hand. This game became popular
    after appearing on the television show 'Big Bang Theory'.
    The rules for the game can be summarised as follows:

    - Scissors cuts Paper
    - Paper covers Rock
    - Rock crushes Lizard
    - Lizard poisons Spock
    - Spock smashes Scissors
    - Scissors decapitates Lizard
    - Lizard eats Paper
    - Paper disproves Spock
    - Spock vaporizes Rock
    - (and as it always has) Rock crushes Scissors

    This can be modeled as a zero sum extensive form game with the following
    tree:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.RPSLSp()
        p = RPSLSp.plot(view_info_sets=True)
        sphinx_plot(p)

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.RPSLSp()
        sage: g
        Rock-Paper-Scissors-Lizard-Spock - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['c', 'd', 'b', 'e', 'f'], 'b': ['Leaf 1', 'Leaf 2', 'Leaf 3', 'Leaf 4', 'Leaf 5'], 'c': ['Leaf 6', 'Leaf 7', 'Leaf 8', 'Leaf 9', 'Leaf 10'], 'd': ['Leaf 11', 'Leaf 12', 'Leaf 13', 'Leaf 14', 'Leaf 15'], 'e': ['Leaf 16', 'Leaf 17', 'Leaf 18', 'Leaf 19', 'Leaf 20'], 'f': ['Leaf 21', 'Leaf 22', 'Leaf 23', 'Leaf 24', 'Leaf 25']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[0, 1], [1, 0], [1, 0], [0, 0], [0, 1], [1, 0], [0, 1], [1, 0], [0, 1], [0, 1],  [0, 0], [1, 0], [1, 0], [0, 1], [1, 0], [1, 0], [0, 1], [0, 0], [0, 0], [1, 0], [0, 1], [0, 0], [0, 1], [1, 0], [0, 1]]
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'L': 0.2, 'P': 0.2, 'R': 0.2, 'S': 0.2, 'Sp': 0.2}), (('c', 'd', 'b', 'e', 'f'), {'L': 0.2, 'P': 0.2, 'R': 0.2, 'S': 0.2, 'Sp': 0.2}))]
        sage: N == d
        True
    """
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_2 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_3 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_4 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_5 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_6 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_7 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_8 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_9 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_10 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_11 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_12 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_13 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_14 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_15 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_16 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_17 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_18 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_19 = EFG_Leaf({player_1: 0, player_2: 0})
    leaf_20 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_21 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_22 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_23 = EFG_Leaf({player_1: 1, player_2: 0})
    leaf_24 = EFG_Leaf({player_1: 0, player_2: 1})
    leaf_25 = EFG_Leaf({player_1: 0, player_2: 0})
    node_6 = EFG_Node(player_2, {'R': leaf_21, 'P': leaf_22, 'S': leaf_23, 'L': leaf_24, 'Sp': leaf_25}, 'f')
    node_5 = EFG_Node(player_2, {'R': leaf_16, 'P': leaf_17, 'S': leaf_18, 'L': leaf_19, 'Sp': leaf_20}, 'e')
    node_4 = EFG_Node(player_2, {'R': leaf_11, 'P': leaf_12, 'S': leaf_13, 'L': leaf_14, 'Sp': leaf_15}, 'd')
    node_3 = EFG_Node(player_2, {'R': leaf_6, 'P': leaf_7, 'S': leaf_8, 'L': leaf_9, 'Sp': leaf_10}, 'c')
    node_2 = EFG_Node(player_2, {'R': leaf_1, 'P': leaf_2, 'S': leaf_3, 'L': leaf_4, 'Sp': leaf_5}, 'b')
    node_1 = EFG_Node(player_1, {'R': node_2, 'P': node_3, 'S': node_4, 'L': node_5, 'Sp': node_6}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3, node_4, node_5, node_6])
    g.rename('Rock-Paper-Scissors-Lizard-Spock - ' + repr(g))
    return g


def Chicken(A=0, a=0, B=1, b=-1, C=-1, c=1, D=-10, d=-10):
    r"""
    Return a Chicken game.

    Consider two drivers locked in a fierce battle for pride. They drive
    towards a cliff and the winner is declared as the last one to swerve.
    If neither player swerves they will both fall off the cliff.

    An example of this can be modeled as an extensive form game tree using the following plot:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.Chicken()
        p = Chicken.plot(view_info_sets=True)
        sphinx_plot(p)

    Where `A < B, D < C` and `a < c, d < b` but with the extra
    condition that `A > C` and `a > b`.

    There are three Nash equilibria:

        1. The second player swerving.
        2. The first player swerving.
        3. Both players swerving with 1 out of 10 times.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.Chicken()
        sage: g
        Chicken - Anti coordination game - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['b', 'c'], 'b': ['Leaf 1', 'Leaf 2'], 'c': ['Leaf 3', 'Leaf 4']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[-10, -10], [1, -1], [-1, 1], [0, 0]]
        sage: N = g.obtain_nash()  # optional - gambit
        sage: d = [((('a',), {'Y': 0.0, 'X': 1.0}), (('b', 'c'), {'Y': 1.0, 'X': 0.0})), ((('a',), {'Y': 0.1, 'X': 0.9}), (('b', 'c'), {'Y': 0.1, 'X': 0.9})), ((('a',), {'Y': 1.0, 'X': 0.0}), (('b', 'c'), {'Y': 0.0, 'X': 1.0}))]
        sage: N == d
        True

    Non default values can be passed::

        sage: g = game_theory.extensive_form_game.Chicken(A=0, a=0, B=2,
        ....:                               b=-1, C=-1, c=2, D=-100, d=-100)
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[-100, -100], [2, -1], [-1, 2], [0, 0]]

    Note that an error is returned if the defining inequalities are not obeyed
    `B > A > C > D` and `c > a > b > d`::

        sage: g = game_theory.extensive_form_game.Chicken(A=8, a=3, B=4, b=2,
        ....:                                               C=2, c=8, D=1, d=0)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a game of chicken must be of the form B > A > C > D and c > a > b > d
    """
    if not (B > A > C > D and c > a > b > d):
        raise TypeError("the input values for a game of chicken must be of the form B > A > C > D and c > a > b > d")
    g = AntiCoordinationGame(A=A, a=a, B=B, b=b, C=C, c=c, D=D, d=d)
    g.rename('Chicken - ' + repr(g))
    return g

def TravellersDilemma():
    r"""
    Return a Travellers dilemma game.

    An airline loses two suitcases belonging to two different travelers. Both
    suitcases happen to be identical and contain identical antiques. An
    airline manager tasked to settle the claims of both travelers explains
    that the airline is liable for a maximum of 10 per suitcase, and in order
    to determine an honest appraised value of the antiques the manager
    separates both travelers so they can't confer, and asks them to write down
    the amount of their value at no less than 2 and no larger than 10. He
    also tells them that if both write down the same number, he will treat
    that number as the true dollar value of both suitcases and reimburse both
    travelers that amount.

    However, if one writes down a smaller number than the other, this smaller
    number will be taken as the true dollar value, and both travelers will
    receive that amount along with a bonus/malus: 2 extra will be paid to the
    traveler who wrote down the lower value and a 2 deduction will be taken
    from the person who wrote down the higher amount. The challenge is: what
    strategy should both travelers follow to decide the value they should
    write down? [Basu]_

    This can be modeled as an extensive form game using the following tree:

    .. PLOT::
        :width: 500 px

        game_theory.extensive_form_game.TravellersDilemma()
        p = TravellersDilemma.plot(view_info_sets=True)
        sphinx_plot(p)


    There is a single Nash equilibrium to this game resulting in
    both players naming the smallest possible value.

    This can be implemented in Sage using the following::

        sage: g = game_theory.extensive_form_game.TravellersDilemma()
        sage: g
        Travellers dilemma - Extensive Form Game with the following underlying tree: {...}
        sage: expected_dictionary = {'a': ['d', 'g', 'e', 'j', 'i', 'h', 'f', 'b', 'c'], 'b': ['Leaf 1', 'Leaf 2', 'Leaf 3', 'Leaf 4', 'Leaf 5', 'Leaf 6', 'Leaf 7', 'Leaf 8', 'Leaf 9'], 'c': ['Leaf 10', 'Leaf 11', 'Leaf 12', 'Leaf 13', 'Leaf 14', 'Leaf 15', 'Leaf 16', 'Leaf 17', 'Leaf 18'], 'd': ['Leaf 19', 'Leaf 20', 'Leaf 21', 'Leaf 22', 'Leaf 23', 'Leaf 24', 'Leaf 25', 'Leaf 26', 'Leaf 27'], 'e': ['Leaf 28', 'Leaf 29', 'Leaf 30', 'Leaf 31', 'Leaf 32', 'Leaf 33', 'Leaf 34', 'Leaf 35', 'Leaf 36'], 'f': ['Leaf 37', 'Leaf 38', 'Leaf 39', 'Leaf 40', 'Leaf 41', 'Leaf 42', 'Leaf 43', 'Leaf 44', 'Leaf 45'], 'g': ['Leaf 46', 'Leaf 47', 'Leaf 48', 'Leaf 49', 'Leaf 50', 'Leaf 51', 'Leaf 52', 'Leaf 53', 'Leaf 54'], 'h': ['Leaf 55', 'Leaf 56', 'Leaf 57', 'Leaf 58', 'Leaf 59', 'Leaf 60', 'Leaf 61', 'Leaf 62', 'Leaf 63'], 'i': ['Leaf 64', 'Leaf 65', 'Leaf 66', 'Leaf 67', 'Leaf 68', 'Leaf 69', 'Leaf 70', 'Leaf 71', 'Leaf 72'], 'j': ['Leaf 73', 'Leaf 74', 'Leaf 75', 'Leaf 76', 'Leaf 77', 'Leaf 78', 'Leaf 79', 'Leaf 80', 'Leaf 81']}
        sage: tree_dictionary = {key.name:[e.name for e in g.tree_dictionary[key]] for key in g.tree_dictionary}
        sage: tree_dictionary == expected_dictionary
        True
        sage: N = g.obtain_nash() # optional - gambit
        sage: d = [((('a',), {'Seven': 0.0, 'Ten': 0.0, 'Six': 0.0, 'Two': 1.0, 'Three': 0.0, 'Four': 0.0, 'Nine': 0.0, 'Eight': 0.0, 'Five': 0.0}), (('d', 'g', 'e', 'i', 'h', 'f', 'b', 'c', 'j'), {'Seven': 0.0, 'Ten': 0.0, 'Six': 0.0, 'Two': 1.0, 'Three': 0.0, 'Four': 0.0, 'Nine': 0.0, 'Eight': 0.0, 'Five': 0.0}))]
        sage: N == d
        True
        sage: utilities = [d.values() for d in [e.payoffs for e in g.leafs]]
        sage: utilities
        [[4, 0], [5, 1], [5, 1], [5, 1], [5, 1], [5, 1], [5, 1], [5, 1], [0, 4], [3, 3], [4, 4], [4, 0], [6, 2], [6, 2], [6, 2], [6, 2], [6, 2], [6, 2], [0, 4], [1, 5], [2, 6], [7, 3], [4, 0], [5, 5], [7, 3], [7, 3], [7, 3], [7, 3], [0, 4], [1, 5], [2, 6], [8, 4], [3, 7], [4, 0], [8, 4], [8, 4], [8, 4], [6, 6], [0, 4], [1, 5], [2, 6], [7, 7], [3, 7], [9, 5], [4, 0], [9, 5], [9, 5], [4, 8], [0, 4], [1, 5], [2, 6], [5, 9], [3, 7], [10, 6], [10, 6], [4, 0], [8, 8], [4, 8], [0, 4], [1, 5], [2, 6], [5, 9], [3, 7], [11, 7], [9, 9], [6, 10], [4, 0], [4, 8], [0, 4], [1, 5], [2, 6], [5, 9], [3, 7], [10, 10], [7, 11], [6, 10], [4, 8], [2, 2], [0, 4], [1, 5], [4, 0]]
    """
    from sage.game_theory.extensive_form_game import EFG_Player
    from sage.game_theory.extensive_form_game import EFG_Node
    from sage.game_theory.extensive_form_game import EFG_Leaf

    player_1 = EFG_Player('Player one')
    player_2 = EFG_Player('Player two')
    leaf_1 = EFG_Leaf({player_1: 2, player_2: 2})
    leaf_2 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_3 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_4 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_5 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_6 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_7 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_8 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_9 = EFG_Leaf({player_1: 4, player_2: 0})
    leaf_10 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_11 = EFG_Leaf({player_1: 3, player_2: 3})
    leaf_12 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_13 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_14 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_15 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_16 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_17 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_18 = EFG_Leaf({player_1: 5, player_2: 1})
    leaf_19 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_20 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_21 = EFG_Leaf({player_1: 4, player_2: 4})
    leaf_22 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_23 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_24 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_25 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_26 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_27 = EFG_Leaf({player_1: 6, player_2: 2})
    leaf_28 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_29 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_30 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_31 = EFG_Leaf({player_1: 5, player_2: 5})
    leaf_32 = EFG_Leaf({player_1: 7, player_2: 3})
    leaf_33 = EFG_Leaf({player_1: 7, player_2: 3})
    leaf_34 = EFG_Leaf({player_1: 7, player_2: 3})
    leaf_35 = EFG_Leaf({player_1: 7, player_2: 3})
    leaf_36 = EFG_Leaf({player_1: 7, player_2: 3})
    leaf_37 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_38 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_39 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_40 = EFG_Leaf({player_1: 3, player_2: 7})
    leaf_41 = EFG_Leaf({player_1: 6, player_2: 6})
    leaf_42 = EFG_Leaf({player_1: 8, player_2: 4})
    leaf_43 = EFG_Leaf({player_1: 8, player_2: 4})
    leaf_44 = EFG_Leaf({player_1: 8, player_2: 4})
    leaf_45 = EFG_Leaf({player_1: 8, player_2: 4})
    leaf_46 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_47 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_48 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_49 = EFG_Leaf({player_1: 3, player_2: 7})
    leaf_50 = EFG_Leaf({player_1: 4, player_2: 8})
    leaf_51 = EFG_Leaf({player_1: 7, player_2: 7})
    leaf_52 = EFG_Leaf({player_1: 9, player_2: 5})
    leaf_53 = EFG_Leaf({player_1: 9, player_2: 5})
    leaf_54 = EFG_Leaf({player_1: 9, player_2: 5})
    leaf_55 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_56 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_57 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_58 = EFG_Leaf({player_1: 3, player_2: 7})
    leaf_59 = EFG_Leaf({player_1: 4, player_2: 8})
    leaf_60 = EFG_Leaf({player_1: 5, player_2: 9})
    leaf_61 = EFG_Leaf({player_1: 8, player_2: 8})
    leaf_62 = EFG_Leaf({player_1: 10, player_2: 6})
    leaf_63 = EFG_Leaf({player_1: 10, player_2: 6})
    leaf_64 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_65 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_66 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_67 = EFG_Leaf({player_1: 3, player_2: 7})
    leaf_68 = EFG_Leaf({player_1: 4, player_2: 8})
    leaf_69 = EFG_Leaf({player_1: 5, player_2: 9})
    leaf_70 = EFG_Leaf({player_1: 6, player_2: 10})
    leaf_71 = EFG_Leaf({player_1: 9, player_2: 9})
    leaf_72 = EFG_Leaf({player_1: 11, player_2: 7})
    leaf_73 = EFG_Leaf({player_1: 0, player_2: 4})
    leaf_74 = EFG_Leaf({player_1: 1, player_2: 5})
    leaf_75 = EFG_Leaf({player_1: 2, player_2: 6})
    leaf_76 = EFG_Leaf({player_1: 3, player_2: 7})
    leaf_77 = EFG_Leaf({player_1: 4, player_2: 8})
    leaf_78 = EFG_Leaf({player_1: 5, player_2: 9})
    leaf_79 = EFG_Leaf({player_1: 6, player_2: 10})
    leaf_80 = EFG_Leaf({player_1: 7, player_2: 11})
    leaf_81 = EFG_Leaf({player_1: 10, player_2: 10})
    node_10 = EFG_Node(player_2, {'Two': leaf_73, 'Three': leaf_74, 'Four': leaf_75, 'Five': leaf_76, 'Six': leaf_77, 'Seven': leaf_78, 'Eight': leaf_79, 'Nine': leaf_80, 'Ten': leaf_81}, 'j')
    node_9 = EFG_Node(player_2, {'Two': leaf_64, 'Three': leaf_65, 'Four': leaf_66, 'Five': leaf_67, 'Six': leaf_68, 'Seven': leaf_69, 'Eight': leaf_70, 'Nine': leaf_71, 'Ten': leaf_72}, 'i')
    node_8 = EFG_Node(player_2, {'Two': leaf_55, 'Three': leaf_56, 'Four': leaf_57, 'Five': leaf_58, 'Six': leaf_59, 'Seven': leaf_60, 'Eight': leaf_61, 'Nine': leaf_62, 'Ten': leaf_63}, 'h')
    node_7 = EFG_Node(player_2, {'Two': leaf_46, 'Three': leaf_47, 'Four': leaf_48, 'Five': leaf_49, 'Six': leaf_50, 'Seven': leaf_51, 'Eight': leaf_52, 'Nine': leaf_53, 'Ten': leaf_54}, 'g')
    node_6 = EFG_Node(player_2, {'Two': leaf_37, 'Three': leaf_38, 'Four': leaf_39, 'Five': leaf_40, 'Six': leaf_41, 'Seven': leaf_42, 'Eight': leaf_43, 'Nine': leaf_44, 'Ten': leaf_45}, 'f')
    node_5 = EFG_Node(player_2, {'Two': leaf_28, 'Three': leaf_29, 'Four': leaf_30, 'Five': leaf_31, 'Six': leaf_32, 'Seven': leaf_33, 'Eight': leaf_34, 'Nine': leaf_35, 'Ten': leaf_36}, 'e')
    node_4 = EFG_Node(player_2, {'Two': leaf_19, 'Three': leaf_20, 'Four': leaf_21, 'Five': leaf_22, 'Six': leaf_23, 'Seven': leaf_24, 'Eight': leaf_25, 'Nine': leaf_26, 'Ten': leaf_27}, 'd')
    node_3 = EFG_Node(player_2, {'Two': leaf_10, 'Three': leaf_11, 'Four': leaf_12, 'Five': leaf_13, 'Six': leaf_14, 'Seven': leaf_15, 'Eight': leaf_16, 'Nine': leaf_17, 'Ten': leaf_18}, 'c')
    node_2 = EFG_Node(player_2, {'Two': leaf_1, 'Three': leaf_2, 'Four': leaf_3, 'Five': leaf_4, 'Six': leaf_5, 'Seven': leaf_6, 'Eight': leaf_7, 'Nine': leaf_8, 'Ten': leaf_9}, 'b')
    node_1 = EFG_Node(player_1, {'Two': node_2, 'Three': node_3, 'Four': node_4, 'Five': node_5, 'Six': node_6, 'Seven': node_7, 'Eight': node_8, 'Nine': node_9, 'Ten': node_10}, 'a')
    g = ExtensiveFormGame(node_1)
    g.set_info_set([node_2, node_3, node_4, node_5, node_6, node_7, node_8, node_9, node_10])
    g.rename('Travellers dilemma - ' + repr(g))
    return g
