r"""
A catalog of normal form games.

This allows us to construct common games directly::

    sage: g = game_theory.normal_form_games.PrisonersDilemma()
    sage: g
    Prisoners dilemma - Normal Form Game with the following utilities: ...

We can then immediately obtain the Nash equilibrium for this game::

    sage: g.obtain_nash()
    [[(0, 1), (0, 1)]]

When we test whether the game is actually the one in question, sometimes we will
build a dictionary to test it, since the printed representation can be
platform-dependent, like so::

    sage: d = {(0, 0): [-2, -2], (0, 1): [-5, 0], (1, 0): [0, -5], (1, 1): [-4, -4]}
    sage: g == d
    True

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

.. [Savani] R. Savani, B. von Stengel.
   *Unit Vector Games*.
   International Journal of Economic Theory (2015).

.. [McLennan] A. McLennan, R. Tourky.
   *Imitation games and computation*
   Games and Economic Behavior (2010)

AUTHOR:

- James Campbell and Vince Knight (06-2014)
"""
from sage.game_theory.normal_form_game import NormalFormGame
from sage.matrix.constructor import matrix


def PrisonersDilemma(R=-2, P=-4, S=-5, T=0):
    r"""
    Return a Prisoners dilemma game.

    Assume two thieves have been caught by the police
    and separated for questioning.
    If both thieves cooperate and do not divulge any information they will
    each get a short sentence.
    If one defects he/she is offered a deal while the other thief will get a
    long sentence.
    If they both defect they both get a medium length sentence.

    This can be modeled as a normal form game using the following two matrices
    [Webb]_:

    .. MATH::

        A = \begin{pmatrix}
            R&S\\
            T&P\\
            \end{pmatrix}


        B = \begin{pmatrix}
            R&T\\
            S&P\\
            \end{pmatrix}

    Where `T > R > P > S`.

    - `R` denotes the reward received for cooperating.
    - `S` denotes the 'sucker' utility.
    - `P` denotes the utility for punishing the other player.
    - `T` denotes the temptation payoff.

    An often used version [Webb]_ is the following:

    .. MATH::

        A = \begin{pmatrix}
            -2&-5\\
            0&-4\\
            \end{pmatrix}


        B = \begin{pmatrix}
            -2&0\\
            -5&-4\\
            \end{pmatrix}

    There is a single Nash equilibrium for this at which both thieves defect.
    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.PrisonersDilemma()
        sage: g
        Prisoners dilemma - Normal Form Game with the following utilities: ...
        sage: d = {(0, 0): [-2, -2], (0, 1): [-5, 0], (1, 0): [0, -5],
        ....:      (1, 1): [-4, -4]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (0, 1)]]

    Note that we can pass other values of R, P, S, T::

        sage: g = game_theory.normal_form_games.PrisonersDilemma(R=-1, P=-2, S=-3, T=0)
        sage: g
        Prisoners dilemma - Normal Form Game with the following utilities:...
        sage: d = {(0, 1): [-3, 0], (1, 0): [0, -3],
        ....:      (0, 0): [-1, -1], (1, 1): [-2, -2]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (0, 1)]]

    If we pass values that fail the defining requirement: `T > R > P > S`
    we get an error message::

        sage: g = game_theory.normal_form_games.PrisonersDilemma(R=-1, P=-2, S=0, T=5)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a Prisoners Dilemma must be
        of the form T > R > P > S

    """
    if not (T > R > P > S):
        raise TypeError("the input values for a Prisoners Dilemma must be of the form T > R > P > S")
    from sage.matrix.constructor import matrix
    A = matrix([[R, S], [T, P]])
    g = NormalFormGame([A, A.transpose()])
    g.rename('Prisoners dilemma - ' + repr(g))
    return g


def CoordinationGame(A=10, a=5, B=0, b=0, C=0, c=0, D=5, d=10):
    r"""
    Return a 2 by 2 Coordination Game.

    A coordination game is a particular type of game where the pure Nash
    equilibrium is for the players to pick the same strategies [Webb]_.

    In general these are represented as a normal form game using the
    following two matrices:

    .. MATH::

        A = \begin{pmatrix}
            A&C\\
            B&D\\
            \end{pmatrix}

        B = \begin{pmatrix}
            a&c\\
            b&d\\
            \end{pmatrix}

    Where `A > B, D > C` and `a > c, d > b`.

    An often used version is the following:

    .. MATH::

        A = \begin{pmatrix}
            10&0\\
            0&5\\
            \end{pmatrix}


        B = \begin{pmatrix}
            5&0\\
            0&10\\
            \end{pmatrix}

    This is the default version of the game created by this function::

        sage: g = game_theory.normal_form_games.CoordinationGame()
        sage: g
        Coordination game - Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [0, 0], (1, 0): [0, 0],
        ....:      (0, 0): [10, 5], (1, 1): [5, 10]}
        sage: g == d
        True

    There are two pure Nash equilibria and one mixed::

        sage: g.obtain_nash()
        [[(0, 1), (0, 1)], [(2/3, 1/3), (1/3, 2/3)], [(1, 0), (1, 0)]]

    We can also pass different values of the input parameters::

        sage: g = game_theory.normal_form_games.CoordinationGame(A=9, a=6,
        ....:                                B=2, b=1, C=0, c=1, D=4, d=11)
        sage: g
        Coordination game - Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [0, 1], (1, 0): [2, 1],
        ....:     (0, 0): [9, 6], (1, 1): [4, 11]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (0, 1)], [(2/3, 1/3), (4/11, 7/11)], [(1, 0), (1, 0)]]

    Note that an error is returned if the defining inequalities are
    not obeyed `A > B, D > C` and `a > c, d > b`::

        sage: g = game_theory.normal_form_games.CoordinationGame(A=9, a=6,
        ....:                               B=0, b=1, C=2, c=10, D=4, d=11)
        Traceback (most recent call last):
        ...
        TypeError: the input values for a Coordination game must
                        be of the form A > B, D > C, a > c and d > b
    """
    if not (A > B  and  D > C and a > c and d > b):
        raise TypeError("the input values for a Coordination game must be of the form A > B, D > C, a > c and d > b")
    from sage.matrix.constructor import matrix
    A = matrix([[A, C], [B, D]])
    B = matrix([[a, c], [b, d]])
    g = NormalFormGame([A, B])
    g.rename('Coordination game - ' + repr(g))
    return g


def BattleOfTheSexes():
    r"""
    Return a Battle of the Sexes game.

    Consider two payers: Amy and Bob.
    Amy prefers to play video games and Bob prefers to
    watch a movie. They both however want to spend their evening
    together.
    This can be modeled as a normal form game using the following two matrices
    [Webb]_:

    .. MATH::

        A = \begin{pmatrix}
            3&1\\
            0&2\\
            \end{pmatrix}


        B = \begin{pmatrix}
            2&1\\
            0&3\\
            \end{pmatrix}

    This is a particular type of Coordination Game.
    There are three Nash equilibria:

    1. Amy and Bob both play video games;
    2. Amy and Bob both watch a movie;
    3. Amy plays video games 75% of the time and Bob watches a movie 75% of the time.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.BattleOfTheSexes()
        sage: g
        Battle of the sexes - Coordination game -
         Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [1, 1], (1, 0): [0, 0], (0, 0): [3, 2], (1, 1): [2, 3]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (0, 1)], [(3/4, 1/4), (1/4, 3/4)], [(1, 0), (1, 0)]]
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
    hare by himself, but a hare is worth less than a stag.

    This can be modeled as a normal form game using the following two matrices
    [Skyrms]_:

    .. MATH::

        A = \begin{pmatrix}
            5&0\\
            4&2\\
            \end{pmatrix}


        B = \begin{pmatrix}
            5&4\\
            0&2\\
            \end{pmatrix}

    This is a particular type of Coordination Game.
    There are three Nash equilibria:

        1. Both friends hunting the stag.
        2. Both friends hunting the hare.
        3. Both friends hunting the stag 2/3rds of the time.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.StagHunt()
        sage: g
        Stag hunt - Coordination game -
         Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [0, 4], (1, 0): [4, 0],
        ....:      (0, 0): [5, 5], (1, 1): [2, 2]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (0, 1)], [(2/3, 1/3), (2/3, 1/3)], [(1, 0), (1, 0)]]

    """
    g = CoordinationGame(A=5, a=5, B=4, b=0, C=0, c=4, D=2, d=2)
    g.rename('Stag hunt - ' + repr(g))
    return g


def AntiCoordinationGame(A=3, a=3, B=5, b=1, C=1, c=5, D=0, d=0):
    r"""
    Return a 2 by 2 AntiCoordination Game.

    An anti coordination game is a particular type of game where the pure Nash
    equilibria is for the players to pick different strategies strategies.

    In general these are represented as a normal form game using the
    following two matrices:

    .. MATH::

        A = \begin{pmatrix}
            A&C\\
            B&D\\
            \end{pmatrix}

        B = \begin{pmatrix}
            a&c\\
            b&d\\
            \end{pmatrix}

    Where `A < B, D < C` and `a < c, d < b`.

    An often used version is the following:

    .. MATH::

        A = \begin{pmatrix}
            3&1\\
            5&0\\
            \end{pmatrix}


        B = \begin{pmatrix}
            3&5\\
            1&0\\
            \end{pmatrix}

    This is the default version of the game created by this function::

        sage: g = game_theory.normal_form_games.AntiCoordinationGame()
        sage: g
        Anti coordination game - Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [1, 5], (1, 0): [5, 1],
        ....:     (0, 0): [3, 3], (1, 1): [0, 0]}
        sage: g == d
        True

    There are two pure Nash equilibria and one mixed::

        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(1/3, 2/3), (1/3, 2/3)], [(1, 0), (0, 1)]]

    We can also pass different values of the input parameters::

        sage: g = game_theory.normal_form_games.AntiCoordinationGame(A=2, a=3,
        ....:                                     B=4, b=2, C=2, c=8, D=1, d=0)
        sage: g
        Anti coordination game - Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [2, 8], (1, 0): [4, 2],
        ....:     (0, 0): [2, 3], (1, 1): [1, 0]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(2/7, 5/7), (1/3, 2/3)], [(1, 0), (0, 1)]]

    Note that an error is returned if the defining inequality is
    not obeyed `A > B, D > C` and `a > c, d > b`::

        sage: g = game_theory.normal_form_games.AntiCoordinationGame(A=8, a=3,
        ....:                                     B=4, b=2, C=2, c=8, D=1, d=0)
        Traceback (most recent call last):
        ...
        TypeError: the input values for an Anti coordination game must be of the form A < B, D < C, a < c and d < b
    """
    if not (A < B  and  D < C and a < c and d < b):
        raise TypeError("the input values for an Anti coordination game must be of the form A < B, D < C, a < c and d < b")
    from sage.matrix.constructor import matrix
    A = matrix([[A, C], [B, D]])
    B = matrix([[a, c], [b, d]])
    g = NormalFormGame([A, B])
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
      `c>v`.

    This can be modeled as a normal form game using the following two matrices
    [Webb]_:

    .. MATH::

        A = \begin{pmatrix}
            v/2-c&v\\
            0&v/2\\
            \end{pmatrix}


        B = \begin{pmatrix}
            v/2-c&0\\
            v&v/2\\
            \end{pmatrix}

    Here are the games with the default values of `v=2` and `c=3`.

    .. MATH::

        A = \begin{pmatrix}
            -2&2\\
            0&1\\
            \end{pmatrix}


        B = \begin{pmatrix}
            -2&0\\
            2&1\\
            \end{pmatrix}

    This is a particular example of an anti coordination game.
    There are three Nash equilibria:

        1. One bird acts like a Hawk and the other like a Dove.
        2. Both birds mix being a Hawk and a Dove

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.HawkDove()
        sage: g
        Hawk-Dove - Anti coordination game -
         Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [2, 0], (1, 0): [0, 2],
        ....:     (0, 0): [-2, -2], (1, 1): [1, 1]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(1/3, 2/3), (1/3, 2/3)], [(1, 0), (0, 1)]]

        sage: g = game_theory.normal_form_games.HawkDove(v=1, c=3)
        sage: g
        Hawk-Dove - Anti coordination game -
         Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [1, 0], (1, 0): [0, 1],
        ....:     (0, 0): [-5/2, -5/2], (1, 1): [1/2, 1/2]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(1/6, 5/6), (1/6, 5/6)], [(1, 0), (0, 1)]]

    Note that an error is returned if the defining inequality is not obeyed
    `c < v`:

        sage: g = game_theory.normal_form_games.HawkDove(v=5, c=1)
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

    This can be modeled as a normal form game using the following two matrices
    [McMillan]_ (we assume that the dominant pig's utilities are given by
    `A`):

    .. MATH::

        A = \begin{pmatrix}
            3&1\\
            6&0\\
            \end{pmatrix}


        B = \begin{pmatrix}
            1&4\\
            -1&0\\
            \end{pmatrix}

    There is a single Nash equilibrium at which the dominant pig pushes the
    lever and the subservient pig does not.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.Pigs()
        sage: g
        Pigs - Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [1, 4], (1, 0): [6, -1],
        ....:     (0, 0): [3, 1], (1, 1): [0, 0]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(1, 0), (0, 1)]]

    """
    from sage.matrix.constructor import matrix
    A = matrix([[3, 1], [6, 0]])
    B = matrix([[1, 4], [-1, 0]])
    g = NormalFormGame([A, B])
    g.rename('Pigs - ' + repr(g))
    return g


def MatchingPennies():
    r"""
    Return a Matching Pennies game.

    Consider two players who can choose to display a coin either Heads
    facing up or Tails facing up.
    If both players show the same face then player 1 wins,
    if not then player 2 wins.

    This can be modeled as a zero sum normal form game with the following
    matrix [Webb]_:

    .. MATH::

        A = \begin{pmatrix}
            1&-1\\
            -1&1\\
            \end{pmatrix}

    There is a single Nash equilibria at which both players randomly
    (with equal probability) pick heads or tails.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.MatchingPennies()
        sage: g
        Matching pennies - Normal Form Game with the following utilities: ...
        sage: d ={(0, 1): [-1, 1], (1, 0): [-1, 1],
        ....:     (0, 0): [1, -1], (1, 1): [1, -1]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(1/2, 1/2), (1/2, 1/2)]]
    """
    from sage.matrix.constructor import matrix
    A = matrix([[1, -1], [-1, 1]])
    g = NormalFormGame([A])
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
    the tie.

    This can be modeled as a zero sum normal form game with the following
    matrix [Webb]_:

    .. MATH::

        A = \begin{pmatrix}
            0 & -1 & 1\\
            1 & 0 & -1\\
            -1 & 1 & 0\\
            \end{pmatrix}

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.RPS()
        sage: g
        Rock-Paper-Scissors - Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [-1, 1], (1, 2): [-1, 1], (0, 0): [0, 0],
        ....:      (2, 1): [1, -1], (1, 1): [0, 0], (2, 0): [-1, 1],
        ....:      (2, 2): [0, 0], (1, 0): [1, -1], (0, 2): [1, -1]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)]]
    """
    from sage.matrix.constructor import matrix
    A = matrix([[0, -1, 1], [1, 0, -1], [-1, 1, 0]])
    g = NormalFormGame([A])
    g.rename('Rock-Paper-Scissors - ' + repr(g))
    return g

def RPSLS():
    r"""
    Return a Rock-Paper-Scissors-Lizard-Spock game.

    `Rock-Paper-Scissors-Lizard-Spock
    <http://www.samkass.com/theories/RPSSL.html>`_ is an extension of
    Rock-Paper-Scissors.
    It is a zero sum game usually played between two
    players where each player simultaneously forms one of three
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

    This can be modeled as a zero sum normal form game with the following
    matrix:

    .. MATH::

        A = \begin{pmatrix}
            0 & -1 & 1 & 1 & -1\\
            1 & 0 & -1 & -1 & 1\\
            -1 & 1 & 0 & 1 & -1\\
            -1 & 1 & -1 & 0 & 1\\
            1 & -1 & 1 & -1 & 0\\
            \end{pmatrix}

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.RPSLS()
        sage: g
        Rock-Paper-Scissors-Lizard-Spock -
         Normal Form Game with the following utilities: ...
        sage: d = {(1, 3): [-1, 1], (3, 0): [-1, 1], (2, 1): [1, -1],
        ....:      (0, 3): [1, -1], (4, 0): [1, -1], (1, 2): [-1, 1],
        ....:      (3, 3): [0, 0], (4, 4): [0, 0], (2, 2): [0, 0],
        ....:      (4, 1): [-1, 1], (1, 1): [0, 0], (3, 2): [-1, 1],
        ....:      (0, 0): [0, 0], (0, 4): [-1, 1], (1, 4): [1, -1],
        ....:      (2, 3): [1, -1], (4, 2): [1, -1], (1, 0): [1, -1],
        ....:      (0, 1): [-1, 1], (3, 1): [1, -1], (2, 4): [-1, 1],
        ....:      (2, 0): [-1, 1], (4, 3): [-1, 1], (3, 4): [1, -1],
        ....:      (0, 2): [1, -1]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(1/5, 1/5, 1/5, 1/5, 1/5), (1/5, 1/5, 1/5, 1/5, 1/5)]]
    """
    from sage.matrix.constructor import matrix
    A = matrix([[0, -1, 1, 1, -1],
                [1, 0, -1, -1, 1],
                [-1, 1, 0, 1 , -1],
                [-1, 1, -1, 0, 1],
                [1, -1, 1, -1, 0]])
    g = NormalFormGame([A])
    g.rename('Rock-Paper-Scissors-Lizard-Spock - ' + repr(g))
    return g


def Chicken(A=0, a=0, B=1, b=-1, C=-1, c=1, D=-10, d=-10):
    r"""
    Return a Chicken game.

    Consider two drivers locked in a fierce battle for pride. They drive
    towards a cliff and the winner is declared as the last one to swerve.
    If neither player swerves they will both fall off the cliff.

    This can be modeled as a particular type of anti coordination game
    using the following two matrices:

    .. MATH::

        A = \begin{pmatrix}
            A&C\\
            B&D\\
            \end{pmatrix}

        B = \begin{pmatrix}
            a&c\\
            b&d\\
            \end{pmatrix}

    Where `A < B, D < C` and `a < c, d < b` but with the extra
    condition that `A > C` and `a > b`.

    Here are the numeric values used by default [Watson]_:

    .. MATH::

        A = \begin{pmatrix}
            0&-1\\
            1&-10\\
            \end{pmatrix}


        B = \begin{pmatrix}
            0&1\\
            -1&-10\\
            \end{pmatrix}

    There are three Nash equilibria:

        1. The second player swerving.
        2. The first player swerving.
        3. Both players swerving with 1 out of 10 times.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.Chicken()
        sage: g
        Chicken - Anti coordination game -
         Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [-1, 1], (1, 0): [1, -1],
        ....:      (0, 0): [0, 0], (1, 1): [-10, -10]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(9/10, 1/10), (9/10, 1/10)], [(1, 0), (0, 1)]]

    Non default values can be passed::

        sage: g = game_theory.normal_form_games.Chicken(A=0, a=0, B=2,
        ....:                               b=-1, C=-1, c=2, D=-100, d=-100)
        sage: g
        Chicken - Anti coordination game -
         Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [-1, 2], (1, 0): [2, -1],
        ....:      (0, 0): [0, 0], (1, 1): [-100, -100]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 1), (1, 0)], [(99/101, 2/101), (99/101, 2/101)],
         [(1, 0), (0, 1)]]

    Note that an error is returned if the defining inequalities are not obeyed
    `B > A > C > D` and `c > a > b > d`::

        sage: g = game_theory.normal_form_games.Chicken(A=8, a=3, B=4, b=2,
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

def TravellersDilemma(max_value=10):
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
    write down?

    This can be modeled as a normal form game using the following two matrices
    [Basu]_:

    .. MATH::

        A = \begin{pmatrix}
            10 & 7  & 6 & 5 & 4 & 3 & 2 & 1 & 0\\
            11 & 9  & 6 & 5 & 4 & 3 & 2 & 1 & 0\\
            10 & 10 & 8 & 5 & 4 & 3 & 2 & 1 & 0\\
            9  & 9  & 9 & 7 & 4 & 3 & 2 & 1 & 0\\
            8  & 8  & 8 & 8 & 6 & 3 & 2 & 1 & 0\\
            7  & 7  & 7 & 7 & 7 & 5 & 2 & 1 & 0\\
            6  & 6  & 6 & 6 & 6 & 6 & 4 & 1 & 0\\
            5  & 5  & 5 & 5 & 5 & 5 & 5 & 3 & 0\\
            4  & 4  & 4 & 4 & 4 & 4 & 4 & 4 & 2\\
            \end{pmatrix}

        B = \begin{pmatrix}
            10 & 11 & 10 & 9 & 8 & 7 & 6 & 5 & 4\\
            7  & 9  & 10 & 9 & 8 & 7 & 6 & 5 & 4\\
            6  & 6  & 8  & 9 & 8 & 7 & 6 & 5 & 4\\
            5  & 5  & 5  & 7 & 8 & 7 & 6 & 5 & 4\\
            4  & 4  & 4  & 4 & 6 & 7 & 6 & 5 & 4\\
            3  & 3  & 3  & 3 & 3 & 5 & 6 & 5 & 4\\
            2  & 2  & 2  & 2 & 2 & 2 & 4 & 5 & 4\\
            1  & 1  & 1  & 1 & 1 & 1 & 1 & 3 & 4\\
            0  & 0  & 0  & 0 & 0 & 0 & 0 & 0 & 2\\
            \end{pmatrix}


    There is a single Nash equilibrium to this game resulting in
    both players naming the smallest possible value.

    This can be implemented in Sage using the following::

        sage: g = game_theory.normal_form_games.TravellersDilemma()
        sage: g
        Travellers dilemma - Normal Form Game with the following utilities: ...
        sage: d = {(7, 3): [5, 1], (4, 7): [1, 5], (1, 3): [5, 9],
        ....:      (4, 8): [0, 4], (3, 0): [9, 5], (2, 8): [0, 4],
        ....:      (8, 0): [4, 0], (7, 8): [0, 4], (5, 4): [7, 3],
        ....:      (0, 7): [1, 5], (5, 6): [2, 6], (2, 6): [2, 6],
        ....:      (1, 6): [2, 6], (5, 1): [7, 3], (3, 7): [1, 5],
        ....:      (0, 3): [5, 9], (8, 5): [4, 0], (2, 5): [3, 7],
        ....:      (5, 8): [0, 4], (4, 0): [8, 4], (1, 2): [6, 10],
        ....:      (7, 4): [5, 1], (6, 4): [6, 2], (3, 3): [7, 7],
        ....:      (2, 0): [10, 6], (8, 1): [4, 0], (7, 6): [5, 1],
        ....:      (4, 4): [6, 6], (6, 3): [6, 2], (1, 5): [3, 7],
        ....:      (8, 8): [2, 2], (7, 2): [5, 1], (3, 6): [2, 6],
        ....:      (2, 2): [8, 8], (7, 7): [3, 3], (5, 7): [1, 5],
        ....:      (5, 3): [7, 3], (4, 1): [8, 4], (1, 1): [9, 9],
        ....:      (2, 7): [1, 5], (3, 2): [9, 5], (0, 0): [10, 10],
        ....:      (6, 6): [4, 4], (5, 0): [7, 3], (7, 1): [5, 1],
        ....:      (4, 5): [3, 7], (0, 4): [4, 8], (5, 5): [5, 5],
        ....:      (1, 4): [4, 8], (6, 0): [6, 2], (7, 5): [5, 1],
        ....:      (2, 3): [5, 9], (2, 1): [10, 6], (8, 7): [4, 0],
        ....:      (6, 8): [0, 4], (4, 2): [8, 4], (1, 0): [11, 7],
        ....:      (0, 8): [0, 4], (6, 5): [6, 2], (3, 5): [3, 7],
        ....:      (0, 1): [7, 11], (8, 3): [4, 0], (7, 0): [5, 1],
        ....:      (4, 6): [2, 6], (6, 7): [1, 5], (8, 6): [4, 0],
        ....:      (5, 2): [7, 3], (6, 1): [6, 2], (3, 1): [9, 5],
        ....:      (8, 2): [4, 0], (2, 4): [4, 8], (3, 8): [0, 4],
        ....:      (0, 6): [2, 6], (1, 8): [0, 4], (6, 2): [6, 2],
        ....:      (4, 3): [8, 4], (1, 7): [1, 5], (0, 5): [3, 7],
        ....:      (3, 4): [4, 8], (0, 2): [6, 10], (8, 4): [4, 0]}
        sage: g == d
        True
        sage: g.obtain_nash() # optional - lrs
        [[(0, 0, 0, 0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 0, 0, 0, 0, 1)]]

    Note that this command can be used to create travellers dilemma for a
    different maximum value of the luggage. Below is an implementation
    with a maximum value of 5::

        sage: g = game_theory.normal_form_games.TravellersDilemma(5)
        sage: g
        Travellers dilemma - Normal Form Game with the following utilities: ...
        sage: d = {(0, 1): [2, 6], (1, 2): [1, 5], (3, 2): [4, 0],
        ....:      (0, 0): [5, 5], (3, 3): [2, 2], (3, 0): [4, 0],
        ....:      (3, 1): [4, 0], (2, 1): [5, 1], (0, 2): [1, 5],
        ....:      (2, 0): [5, 1], (1, 3): [0, 4], (2, 3): [0, 4],
        ....:      (2, 2): [3, 3], (1, 0): [6, 2], (0, 3): [0, 4],
        ....:      (1, 1): [4, 4]}
        sage: g == d
        True
        sage: g.obtain_nash()
        [[(0, 0, 0, 1), (0, 0, 0, 1)]]

    """
    from sage.matrix.constructor import matrix
    from sage.functions.generalized import sign
    A = matrix([[min(i, j) + 2 * sign(j - i) for j in range(max_value, 1, -1)]
                                             for i in range(max_value, 1, -1)])
    g = NormalFormGame([A, A.transpose()])
    g.rename('Travellers dilemma - ' + repr(g))
    return g

def RandomGame(ring, number_of_strategies, min_bound=None, max_bound=None):
    r"""
    Returns a random N player normal form game with payoffs matrices in the
    specified ring and with number of strategies specified by the
    ``number_of_strategies`` tuple.

    INPUT:
    -  ``ring`` - base ring for entries of the matrix
    -  ``number_of_strategies`` - Tuple; number of strategies available to each
       player.
    -  ``min_bound`` - This corresponds to ``x`` for ``ZZ``, ``num_bound`` for ``QQ`` and ``min``
       for ``RR``.
    -  ``max_bound`` - This corresponds to ``y`` for ``ZZ``, ``den_bound`` for ``QQ`` and ``max``
       for ``RR``.

    .. NOTE ::
        The meaning of `min_bound` and `max_bound` changes depending on the ring being used as this
        function relies on the `random_element` function which is implemented by each ring

    EXAMPLES:
    Here is a basic random 5 by 5 bimatrix game::

        sage: g = game_theory.normal_form_games.RandomGame(ZZ, (5, 5))
        sage: g.payoff_matrices()
        (
        [ -8   0   1   2 -95]  [  2   0  -1   1  -1]
        [ -2   0   1   1  -2]  [-12   0  -1  -1  -1]
        [  4  -6   0  -2   1]  [ -4   5   0   0  -4]
        [ -6  -1   1   1  -3]  [  1   1  -1  -1   1]
        [  1   0   2   0   1], [  0  -3  -2  -2   0]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 0, 1), (0, 0, 0, 0, 1)],
         [(0, 0, 1/4, 0, 3/4), (0, 0, 0, 0, 1)],
         [(0, 0, 1/4, 0, 3/4), (2/3, 1/3, 0, 0, 0)],
         [(0, 3/4, 0, 0, 1/4), (0, 1, 0, 0, 0)],
         [(0, 1, 0, 0, 0), (0, 1, 0, 0, 0)],
         [(1/4, 1/2, 0, 0, 1/4), (0, 1, 0, 0, 0)],
         [(1/2, 1/2, 0, 0, 0), (0, 1, 0, 0, 0)],
         [(11/12, 1/12, 0, 0, 0), (1/7, 0, 0, 6/7, 0)]]

    We can create a random game with differing sizes of strategies::

        sage: g = game_theory.normal_form_games.RandomGame(ZZ, (2, 6))
        sage: g.payoff_matrices()
        (
        [-1  0  1  4  1 14]  [ 1  0 -1 -1 -1  1]
        [-5 -1  2  1 -2  0], [ 4  0  4  1 -1  4]
        )
        sage: g.obtain_nash()
        [[(0, 1), (0, 0, 14/15, 0, 0, 1/15)],
         [(0, 1), (0, 0, 1, 0, 0, 0)],
         [(0, 1), (1/5, 0, 4/5, 0, 0, 0)],
         [(1, 0), (0, 0, 0, 0, 0, 1)],
         [(1, 0), (1, 0, 0, 0, 0, 0)]]

    It is possible to create games that have utilities that are not in ZZ::

        sage: g = game_theory.normal_form_games.RandomGame(QQ, (3, 3))
        sage: g.payoff_matrices()
        (
        [     -3     1/9     3/2]  [   5/2      1     -1]
        [   -1/2    -2/7 -3/1955]  [    -1    1/2     22]
        [      2   -1/14       1], [     0 -1/357      2]
        )
        sage: g.obtain_nash()
        [[(4/11, 0, 7/11), (1/11, 0, 10/11)]]

    We can easily create a random game with more than 2 players::

        sage: g = game_theory.normal_form_games.RandomGame(ZZ, (4, 2, 3))
        sage: g
        Normal Form Game with the following utilities: ...
        sage: expected_utilities = {(0, 0, 0): [0, -1, -2],
        ....:                       (0, 0, 1): [0, -1, -2],
        ....:                       (0, 0, 2): [0, 0, 27],
        ....:                       (0, 1, 0): [-1, 1, 1],
        ....:                       (0, 1, 1): [0, 2, -1],
        ....:                       (0, 1, 2): [1, -1, -2],
        ....:                       (1, 0, 0): [-1, 3, 2],
        ....:                       (1, 0, 1): [-2, -1, -4],
        ....:                       (1, 0, 2): [2, 1, 2],
        ....:                       (1, 1, 0): [1, -2, 0],
        ....:                       (1, 1, 1): [0, 1, -2],
        ....:                       (1, 1, 2): [-2, -2, 0],
        ....:                       (2, 0, 0): [2, 0, 1],
        ....:                       (2, 0, 1): [1, 1, 1],
        ....:                       (2, 0, 2): [-1, -2, 0],
        ....:                       (2, 1, 0): [-1, 2, -1],
        ....:                       (2, 1, 1): [1, -27, 2],
        ....:                       (2, 1, 2): [5, 1, -10],
        ....:                       (3, 0, 0): [2, -1, -1],
        ....:                       (3, 0, 1): [-2, 2, -4],
        ....:                       (3, 0, 2): [0, 0, 0],
        ....:                       (3, 1, 0): [-1, -1, 1],
        ....:                       (3, 1, 1): [-4, 0, 3],
        ....:                       (3, 1, 2): [2, 1, 0]}
        sage: g.utilities  == expected_utilities  # Testing the equality of utilities as ordering of dictionaries is system dependent.
        True
        sage: g = game_theory.normal_form_games.RandomGame(ZZ, (1, 1, 1, 1, 1, 1))
        sage: g
        Normal Form Game with the following utilities: {(0, 0, 0, 0, 0, 0): [-1, -13, 1, 2, 3, 0]}

    The ``max_bound`` and ``min_bound`` variables can be used to control the values of the
    payoff entries in the game::

        sage: g = game_theory.normal_form_games.RandomGame(ZZ, (3, 3), 6, 9)
        sage: g.payoff_matrices()
        (
        [6 6 6]  [7 6 6]
        [6 8 8]  [8 7 8]
        [8 8 8], [6 6 7]
        )
        sage: g = game_theory.normal_form_games.RandomGame(QQ, (3, 3), 6, 9)
        sage: g.payoff_matrices()
        (
        [-2/3  2/5    1]  [ 3/2  1/4 -5/3]
        [   2  2/3 -1/2]  [-2/9   -3  1/4]
        [ 1/9 -1/2 -1/2], [ 5/4 -5/9    0]
        )
        sage: g = game_theory.normal_form_games.RandomGame(RR, (3, 3), 0, 1)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [0.0480  0.981  0.164]
        [ 0.338  0.163  0.673]
        [ 0.257  0.374  0.569]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ 0.461  0.950  0.783]
        [ 0.872 0.0375  0.415]
        [ 0.447  0.448  0.331]
    """
    from sage.combinat.cartesian_product import CartesianProduct
    from sage.rings.all import RR
    if ring == RR :
        if min_bound is None:
            min_bound = -1
        if max_bound is None:
            max_bound = 1
    g = NormalFormGame()
    for number in number_of_strategies:
        g.add_player(number)
    for profile in CartesianProduct(*[range(number) for number in
        number_of_strategies]):
        g[tuple(profile)] = [ring.random_element(min_bound, max_bound) for player in number_of_strategies]
    return g

def RandomZeroSum(ring, m, n,  min_bound = None, max_bound = None):
    r"""
    Creates a random 2 player zero sum game (A, -A) of size `m \times n`, where the payoffs
    of the row players are chosen uniformly at random using the ``min_bound`` and ``max_bound``
    variables.

    INPUT:

    -  ``ring`` - base ring for entries in the matrix
    -  ``m`` - Number of strategies for the row player
    -  ``n`` - Number of strategies for the column player
    -  ``min_bound`` - This corresponds to ``x`` for ``ZZ``, ``num_bound`` for ``QQ`` and ``min``
       for ``RR``.
    -  ``max_bound`` - This corresponds to ``y`` for ``ZZ``, ``den_bound`` for ``QQ`` and ``max``
       for ``RR``.

    .. NOTE ::
        The meaning of `min_bound` and `max_bound` changes depending on the ring being used as this
        function relies on the `random_element` function which is implemented by each ring

    EXAMPLE:
    Here is a random zero sum game with 5 actions per player::

        sage: g = game_theory.normal_form_games.RandomZeroSum(ZZ, 5, 5)
        sage: g.payoff_matrices()
        (
        [ -8   2   0   0   1]  [ 8 -2  0  0 -1]
        [ -1   2   1 -95  -1]  [ 1 -2 -1 95  1]
        [ -2 -12   0   0   1]  [ 2 12  0  0 -1]
        [ -1   1  -1  -2  -1]  [ 1 -1  1  2  1]
        [  4  -4  -6   5   0], [-4  4  6 -5  0]
        )
        sage: g.obtain_nash()
        [[(18/371, 0, 1/7, 275/371, 25/371), (17/106, 51/742, 181/371, 15/53, 0)]]

    We can create a random zero sum game with differing number of strategies for both players::

        sage: g = game_theory.normal_form_games.RandomZeroSum(ZZ, 2, 6)
        sage: g.payoff_matrices()
        (
        [ 0 -2  0  1 -4 -6]  [  0   2   0  -1   4   6]
        [ 1 -1  1  1 -1  1]  [ -1   1  -1  -1   1  -1]
        [-1 -3  1  1  0  0]  [  1   3  -1  -1   0   0]
        [-3  2 -2  0 -2  1]  [  3  -2   2   0   2  -1]
        [ 0 -1  1  0  0  1]  [  0   1  -1   0   0  -1]
        [-1  4 -1  1 -1 14], [  1  -4   1  -1   1 -14]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 0, 5/6, 1/6), (0, 1/6, 0, 0, 5/6, 0)],
         [(0, 0, 0, 0, 5/6, 1/6), (5/12, 1/6, 0, 0, 5/12, 0)]]

    It is also possible to create games where the utilities are not in ZZ::

        sage: g = game_theory.normal_form_games.RandomZeroSum(RR, 5, 5)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ 0.369  0.821  0.190 -0.808 -0.800]
        [ 0.532 -0.852 -0.444 -0.916  0.706]
        [ 0.316 -0.205 -0.214  0.669  0.897]
        [ 0.108 -0.612 -0.308  0.250  0.575]
        [ 0.869 -0.356 -0.552  0.706 -0.324]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [-0.369 -0.821 -0.190  0.808  0.800]
        [-0.532  0.852  0.444  0.916 -0.706]
        [-0.316  0.205  0.214 -0.669 -0.897]
        [-0.108  0.612  0.308 -0.250 -0.575]
        [-0.869  0.356  0.552 -0.706  0.324]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.469629, 0.0, 0.530371, 0.0, 0.0], [0.0, 0.0, 0.78522, 0.21478, 0.0]]]
        sage: g = game_theory.normal_form_games.RandomZeroSum(QQ, 5, 5)
        sage: g.payoff_matrices()
        (
        [    2     0     2 -2/27    -1]  [  -2    0   -2 2/27    1]
        [  1/2    -1   1/2  -1/3    -1]  [-1/2    1 -1/2  1/3    1]
        [  1/4     2     2    -2     1]  [-1/4   -2   -2    2   -1]
        [   -1     0     1    -1     2]  [   1    0   -1    1   -2]
        [   -2 -1/27   2/5 -1/10    -2], [   2 1/27 -2/5 1/10    2]
        )
        sage: g.obtain_nash()
        [[(81/106, 0, 0, 25/106, 0), (0, 0, 0, 81/106, 25/106)]]

    The range of values of the payoff entries can be controlled using the ``min_bound`` and
    ``max_bound`` variables respectively::

        sage: g = game_theory.normal_form_games.RandomZeroSum(ZZ, 5, 5, -2, 6)
        sage: g.payoff_matrices()
        (
        [-2  1  4  1  4]  [ 2 -1 -4 -1 -4]
        [ 4  2  2  0 -1]  [-4 -2 -2  0  1]
        [ 4 -2  3  2  5]  [-4  2 -3 -2 -5]
        [ 2 -2 -2  1  0]  [-2  2  2 -1  0]
        [ 0  4  3  0  0], [ 0 -4 -3  0  0]
        )
        sage: g.obtain_nash()
        [[(0, 0, 1/2, 0, 1/2), (0, 1/4, 0, 3/4, 0)],
         [(1/4, 0, 3/8, 0, 3/8), (0, 1/4, 0, 3/4, 0)]]
        sage: g = game_theory.normal_form_games.RandomZeroSum(RR, 5, 5, -1, 1)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.674  0.0914  -0.380  -0.694 -0.0509]
        [  0.229   0.602  -0.238  -0.727  -0.310]
        [  0.632   0.660   0.146  -0.196  -0.279]
        [  0.556  -0.642   0.667  -0.512  -0.668]
        [  0.151   0.174   0.850  -0.720   0.371]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.674 -0.0914   0.380   0.694  0.0509]
        [ -0.229  -0.602   0.238   0.727   0.310]
        [ -0.632  -0.660  -0.146   0.196   0.279]
        [ -0.556   0.642  -0.667   0.512   0.668]
        [ -0.151  -0.174  -0.850   0.720  -0.371]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.929997, 0.0, 0.070003], [0.0, 0.0, 0.0, 0.553706, 0.446294]]]
        sage: g = game_theory.normal_form_games.RandomZeroSum(QQ, 5, 5, 5, 8)
        sage: g.payoff_matrices()
        (
        [ 3/5  3/5  5/4 -5/8 -4/3]  [-3/5 -3/5 -5/4  5/8  4/3]
        [ 2/7  3/8  2/7 -5/4    0]  [-2/7 -3/8 -2/7  5/4    0]
        [  -1  1/2  1/3  4/5 -3/5]  [   1 -1/2 -1/3 -4/5  3/5]
        [-1/6    1  4/5 -3/2    2]  [ 1/6   -1 -4/5  3/2   -2]
        [  -5   -1 -1/4  1/2    0], [   5    1  1/4 -1/2    0]
        )
        sage: g.obtain_nash()
        [[(15960/45559, 0, 18835/45559, 10764/45559, 0),
          (19410/45559, 0, 0, 17176/45559, 8973/45559)]]

    """
    from sage.rings.all import RR
    if ring == RR :
        if min_bound is None:
            min_bound = -1
        if max_bound is None:
            max_bound = 1
    A = matrix(ring, n, n)
    for i in range(n):
        for j in range(n):
            A[i, j] = ring.random_element(min_bound, max_bound)
    return NormalFormGame([A])

def RandomCovariantGame(m, n, covariance):
    r"""
    Creates a random 2 player covariant game of size `m \times n` where the payoff
    matrices are chosen from a multivariate normal distribution such that the covariance between the
    payoffs of each player is `covariance`.

    INPUT:

    - ``m`` - The number of strategies for the row player.
    - ``m`` - The number of strategies for the column player.
    - ``covariance`` - The covariance between the two payoff matrices.

    .. NOTE ::
        As the covariance value gets close to `1` it is similar to the game `(A, A)`, as it gets
        closer to zero, the game is similar to the payoffs being chosen independently from the
        normal distribution. Finally, when the covariance gets closer to `-1`, it gets closer to
        being a zero-sum game `(A, -A)`.

    EXAMPLE:
    Here is an example of a 5 by 5 bimatrix game::

        sage: import numpy
        sage: numpy.random.seed(1)
        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, -0.5)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  -1.71 -0.0791   -1.90   -1.89  -0.401]
        [  -2.30  0.0872   -1.53  -0.290   0.255]
        [   1.53  -0.530   -1.12  -0.361   0.497]
        [  0.401   0.173   0.575    1.08   -1.07]
        [ -0.278    1.49  -0.362   0.885   0.205]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  1.10 -0.994 -0.401   1.13  0.152]
        [ 0.236 -0.471  0.432 -0.588  0.328]
        [-0.381   1.03  0.438 -0.574 0.0332]
        [-0.797  -1.02 -0.588 -0.850   1.81]
        [-0.610  0.199 -0.275   1.22  0.413]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.0, 0.231912, 0.768088], [0.0, 0.0, 0.0, 0.864012, 0.135988]]]

    As the ``covariance`` parameter gets closer to `-1`, the games which are generated get close to
    being zero-sum::

        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, -0.9)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.371    1.04   0.335  -0.610 -0.0804]
        [   1.02  -0.567  -0.493  -0.763   -2.44]
        [   1.29  0.0399  -0.760   0.484 -0.0539]
        [  0.172 -0.0901  -0.167   0.738   0.134]
        [  -1.13   0.223  -0.395   0.345   0.760]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.214   -1.19 -0.0724    1.03   0.476]
        [ -0.455   0.433   0.459    1.44    1.82]
        [  -1.52   0.352  -0.145  -0.113   0.395]
        [ -0.262   0.274   0.220  -0.569   0.371]
        [   1.21  -0.509   0.430  -0.325  -0.448]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.23597, 0.0, 0.0, 0.52594, 0.23809],
          [0.159792, 0.0, 0.372795, 0.0, 0.467413]]]

    with games generated with ``covariance`` of `-1.0` being zero-sum games::

        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, -1.0)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.447  -0.404    1.09  -0.741   0.266]
        [   1.37  -0.846  -0.351  0.0387   -1.12]
        [ 0.0246   -1.27    1.86   -1.63    1.20]
        [  0.181    1.23  -0.793  -0.521  -0.802]
        [  0.187  -0.869  -0.529 -0.0778  -0.232]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.447   0.404   -1.09   0.741  -0.266]
        [  -1.37   0.846   0.351 -0.0387    1.12]
        [-0.0246    1.27   -1.86    1.63   -1.20]
        [ -0.181   -1.23   0.793   0.521   0.802]
        [ -0.187   0.869   0.529  0.0778   0.232]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.178578, 0.023132, 0.0, 0.232827, 0.565464],
          [0.0, 0.169005, 0.15081, 0.570704, 0.109481]]]

    And as the ``covariance`` gets closer to `1`, the games get closer to coordination games `(A,
    A)`::

        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, 0.9)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.847   -1.50  -0.408   0.160   0.234]
        [  0.727  -0.176  -0.444   -1.01   -1.71]
        [  0.429   -2.41 0.00804   -1.23  -0.591]
        [   1.17  -0.780   0.246  -0.602  -0.257]
        [  -1.21   -2.10  -0.278  -0.681   -1.02]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.242  -0.523  -0.453   0.106  -0.268]
        [  0.281  -0.309  -0.522  -0.914   -2.56]
        [  0.832   -2.52 -0.0932   -1.36  -0.734]
        [   1.31  -0.201 -0.0304  -0.495  -0.290]
        [  -1.05   -1.61   -1.01  -0.493  -0.557]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0]],
         [[0.838161, 0.0, 0.0, 0.161839, 0.0], [0.701958, 0.0, 0.0, 0.298042, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0]]]


    with games generated with ``covariance`` of `1.0` being coordination games::

        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, 1.0)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.401   0.562    1.33    1.65    1.12]
        [  0.327   -1.11    1.24  -0.623   -1.41]
        [  -1.62   -1.56    1.22   0.546   0.700]
        [ -0.243  -0.661   0.120    1.18    1.67]
        [  0.498 0.00189   0.861  -0.619   -1.81]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.401   0.562    1.33    1.65    1.12]
        [  0.327   -1.11    1.24  -0.623   -1.41]
        [  -1.62   -1.56    1.22   0.546   0.700]
        [ -0.243  -0.661   0.120    1.18    1.67]
        [  0.498 0.00189   0.861  -0.619   -1.81]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.0]],
         [[0.480592, 0.0, 0.0, 0.519408, 0.0], [0.0, 0.0, 0.0, 0.542344, 0.457656]],
         [[1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0]]]

    Finally, a ``covariance`` of 0 results in a game where the matrices are drawn independently from
    the normal distribution::

        sage: g = game_theory.normal_form_games.RandomCovariantGame(5, 5, 0)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.345   -2.79   0.366    2.05   0.430]
        [  0.106   0.795   0.134   0.285   0.276]
        [  0.836   0.759  -0.877   -1.44  -0.254]
        [ -0.782  0.0954  0.0607  0.0165   -1.12]
        [ -0.187   0.492 -0.0845   0.417  -0.955]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ -0.231    1.94   -1.04   0.586  -0.607]
        [  -1.53  -0.374    1.20   0.262  -0.733]
        [   1.54   0.885  -0.868    1.23    1.40]
        [ -0.438   0.921   0.211   0.177  0.0809]
        [-0.0568  -0.681  -0.297   0.785   0.586]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0]],
         [[0.21041, 0.035178, 0.754413, 0.0, 0.0],
          [0.665425, 0.054301, 0.0, 0.280274, 0.0]],
         [[0.345828, 0.654172, 0.0, 0.0, 0.0], [0.0, 0.060801, 0.939199, 0.0, 0.0]]]

    """
    from sage.rings.all import RR
    from scipy.stats import multivariate_normal
    cov = [[1.0, covariance], [covariance, 1.0]]
    d = multivariate_normal.rvs(cov = cov, size = int(m)*int(n)).T
    A = matrix(RR, m, n, list(d[0]))
    B = matrix(RR, m, n, list(d[1]))
    return NormalFormGame([A, B])

def RandomUnitVectorGame(ring, n, min_bound=None, max_bound=None):
    r"""
    Returns a random unit vector game, where the payoffs of the column player are chosen
    uniformly at random, and each row of the row players payoff matrix is a unit column vector
    [Savani]_.

    INPUT:

    -  ``ring`` - the base ring of the elements in the payoff matrix
    -  ``n`` - the number of payoff entries for both players
    -  ``min_bound`` - This corresponds to ``x`` for ``ZZ``, ``num_bound`` for ``QQ`` and ``min``
       for ``RR``.
    -  ``max_bound`` - This corresponds to ``y`` for ``ZZ``, ``den_bound`` for ``QQ`` and ``max``
       for ``RR``.

    .. NOTE ::
        The meaning of `min_bound` and `max_bound` changes depending on the ring being used as this
        function relies on the `random_element` function which is implemented by each ring

    EXAMPLES:
    Here is a basic random 5 by 5 unit vector bimatrix game::

        sage: g = game_theory.normal_form_games.RandomUnitVectorGame(ZZ, 5)
        sage: g.payoff_matrices()
        (
        [0 0 0 1 1]  [ -8   2   0   0   1]
        [1 0 0 0 0]  [ -1   2   1 -95  -1]
        [0 0 0 0 0]  [ -2 -12   0   0   1]
        [0 1 0 0 0]  [ -1   1  -1  -2  -1]
        [0 0 1 0 0], [  4  -4  -6   5   0]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 1, 0), (0, 1, 0, 0, 0)]]

    The payoff entries in the game can also be in other rings::

        sage: g = game_theory.normal_form_games.RandomUnitVectorGame(QQ, 5)
        sage: g.payoff_matrices()
        (
        [0 0 0 1 1]  [      0      -1       0      -1       1]
        [0 0 0 0 0]  [   -1/4      -1   -1/14    -1/5      -4]
        [1 0 0 0 0]  [      0       4    -1/2    -1/4      -3]
        [0 1 1 0 0]  [    5/2     1/9       1     3/2      -1]
        [0 0 0 0 0], [   -1/2      -1    -2/7     1/2 -3/1955]
        )
        sage: g.obtain_nash()
        [[(0, 0, 43/115, 72/115, 0), (1/2, 1/2, 0, 0, 0)],
         [(4528903/11961647, 320790/11961647, 1899924/11961647, 0, 5212030/11961647),
          (0, 0, 0, 0, -1)],
         [(4528903/11961647, 320790/11961647, 1899924/11961647, 0, 5212030/11961647),
          (0, 1, 0, 0, -1)],
         [(381/505, 0, 106/505, 18/505, 0), (1/3, 1/3, 0, 0, 1/3)],
         [(1, 0, 0, 0, 0), (0, 0, 0, 0, 1)]]
        sage: g = game_theory.normal_form_games.RandomUnitVectorGame(RR, 5)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [0.000 0.000 0.000 0.000 0.000]
        [0.000 0.000 0.000  1.00 0.000]
        [0.000 0.000  1.00 0.000 0.000]
        [ 1.00 0.000 0.000 0.000  1.00]
        [0.000  1.00 0.000 0.000 0.000]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.415  -0.649   0.793   0.451 -0.0583]
        [ -0.254   0.223  -0.358   0.598   0.261]
        [  0.437  -0.383   0.782   0.540  -0.748]
        [ -0.677   0.766   0.780   0.550   0.848]
        [  0.854   0.116   0.378  -0.949   0.587]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.0]],
         [[0.0, 0.0, 0.042447, 0.957553, 0.0], [0.0, 0.0, 0.5, 0.0, 0.5]],
         [[0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0]],
         [[0.0, 0.195053, 0.109765, 0.695182, 0.0],
          [0.0, 0.0, 0.333333, 0.333333, 0.333333]],
         [[0.0, 0.202299, 0.797701, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5, 0.0]],
         [[0.0, 0.469638, 0.0, 0.530362, 0.0], [0.0, 0.0, 0.0, 0.5, 0.5]],
         [[0.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 0.0]]]

    The ranges of the payoff entries in the game can be controlled by using the ``max_bound`` and
    ``min_bound`` parameters::

        sage: g = game_theory.normal_form_games.RandomUnitVectorGame(ZZ, 5, -3, 5)
        sage: g.payoff_matrices()
        (
        [0 1 0 0 1]  [-1  2  0  0  4]
        [0 0 1 1 0]  [ 2 -2  4  0  1]
        [0 0 0 0 0]  [-3 -3  0  1 -3]
        [1 0 0 0 0]  [ 4  1  4 -2 -3]
        [0 0 0 0 0], [ 2 -1  1  3  0]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 1, 0), (1/2, 0, 1/2, 0, 0)],
         [(0, 0, 0, 1, 0), (1, 0, 0, 0, 0)],
         [(0, 1, 0, 0, 0), (0, 0, 1, 0, 0)],
         [(3/7, 4/7, 0, 0, 0), (0, 0, 1/2, 0, 1/2)],
         [(1, 0, 0, 0, 0), (0, 0, 0, 0, 1)]]
        sage: g = game_theory.normal_form_games.RandomUnitVectorGame(RR, 5, -3, 5)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [0.000 0.000  1.00  1.00 0.000]
        [0.000  1.00 0.000 0.000 0.000]
        [ 1.00 0.000 0.000 0.000 0.000]
        [0.000 0.000 0.000 0.000  1.00]
        [0.000 0.000 0.000 0.000 0.000]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ 2.57  4.77 -1.56  2.00 0.551]
        [ 4.11 -2.60  3.14 -2.42 -2.54]
        [-1.77  1.77  1.13 -1.62 -2.03]
        [-2.12  3.50 -2.92  3.30  4.20]
        [ 2.30 0.209 -1.50  4.92 0.111]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.0, 0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0, 1.0]],
         [[0.0, 0.44014, 0.092534, 0.467325, 0.0],
          [0.333333, 0.333333, 0.0, 0.0, 0.333333]],
         [[0.239329, 0.314287, 0.446383, 0.0, 0.0],
          [0.333333, 0.333333, 0.333333, 0.0, 0.0]]]

    """
    from random import randint
    from sage.rings.all import RR, ZZ

    if ring == RR :
        if min_bound is None:
            min_bound = -1
        if max_bound is None:
            max_bound = 1

    B = matrix(ring, n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = ring.random_element(min_bound, max_bound)

    A = matrix(ring, n, n)
    for i in range(n):
        j = ZZ.random_element(0, n)
        A[j, i] = 1

    return NormalFormGame([A, B])

def RandomImitationGame(ring, n, min_bound=None, max_bound=None):
    r"""
    Returns a random imitation game, where the payoffs of the column player are chosen
    uniformly at random and the payoff matrix of the row player is the identity matrix [Savani,
    McLennan]_.

    INPUT:

    -  ``ring`` - the base ring for entries of the matrix
    -  ``n`` - the number of strategies for each player
    -  ``min_bound`` - This corresponds to ``x`` for ``ZZ``, ``num_bound`` for ``QQ`` and ``min``
       for ``RR``.
    -  ``max_bound`` - This corresponds to ``y`` for ``ZZ``, ``den_bound`` for ``QQ`` and ``max``
       for ``RR``.

    .. NOTE ::
        The meaning of `min_bound` and `max_bound` changes depending on the ring being used as this
        function relies on the `random_element` function which is implemented by each ring

    EXAMPLES:
    Here is a random 5 by 5 imitation game::

        sage: g = game_theory.normal_form_games.RandomImitationGame(ZZ, 5)
        sage: g.payoff_matrices()
        (
        [1 0 0 0 0]  [ -8   2   0   0   1]
        [0 1 0 0 0]  [ -1   2   1 -95  -1]
        [0 0 1 0 0]  [ -2 -12   0   0   1]
        [0 0 0 1 0]  [ -1   1  -1  -2  -1]
        [0 0 0 0 1], [  4  -4  -6   5   0]
        )
        sage: g.obtain_nash()
        [[(0, 12/13, 1/13, 0, 0), (0, 1/2, 1/2, 0, 0)],
         [(0, 1, 0, 0, 0), (0, 1, 0, 0, 0)],
         [(31/573, 376/573, 55/573, 0, 37/191), (1/4, 1/4, 1/4, 0, 1/4)],
         [(37/267, 5/534, 0, 577/1068, 111/356), (1/4, 1/4, 0, 1/4, 1/4)],
         [(12/71, 32/71, 0, 0, 27/71), (1/3, 1/3, 0, 0, 1/3)]]

    The payoff entries in the game can also be in other rings::

        sage: g = game_theory.normal_form_games.RandomImitationGame(QQ, 5)
        sage: g.payoff_matrices()
        (
        [1 0 0 0 0]  [    0     0   2/3    -1     1]
        [0 1 0 0 0]  [   -1   1/3     1     0    -1]
        [0 0 1 0 0]  [    0    -1     1  -1/4    -1]
        [0 0 0 1 0]  [-1/14  -1/5    -4     0     4]
        [0 0 0 0 1], [ -1/2  -1/4    -3   5/2   1/9]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 43/115, 72/115), (0, 0, 0, 1/2, 1/2)],
         [(0, 0, 1, 0, 0), (0, 0, 1, 0, 0)],
         [(0, 11012/17913, 9656/53739, 10435/53739, 68/5971), (0, 1/4, 1/4, 1/4, 1/4)],
         [(0, 6198/8725, 0, 383/1745, 612/8725), (0, 1/3, 0, 1/3, 1/3)],
         [(116/1437, 0, 9029/12933, 3605/25866, 235/2874), (1/4, 0, 1/4, 1/4, 1/4)]]
        sage: g = game_theory.normal_form_games.RandomImitationGame(RR, 5)
        sage: A, B = g.payoff_matrices()
        sage: print A.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [ 1.00 0.000 0.000 0.000 0.000]
        [0.000  1.00 0.000 0.000 0.000]
        [0.000 0.000  1.00 0.000 0.000]
        [0.000 0.000 0.000  1.00 0.000]
        [0.000 0.000 0.000 0.000  1.00]
        sage: print B.str(rep_mapping=lambda x : str(x.n(digits=3)))
        [  0.316  -0.205  -0.214   0.669   0.897]
        [  0.108  -0.612  -0.308   0.250   0.575]
        [  0.869  -0.356  -0.552   0.706  -0.324]
        [  0.415  -0.649   0.793   0.451 -0.0583]
        [ -0.254   0.223  -0.358   0.598   0.261]
        sage: eqs = g.obtain_nash()
        sage: [[[round(float(p), 6) for p in st] for st in eq] for eq in eqs]
        [[[0.014305, 0.0, 0.202798, 0.782897, 0.0],
          [0.333333, 0.0, 0.333333, 0.333333, 0.0]]]

    The ranges of the payoff entries in the game can be controlled by using the ``max_bound`` and
    ``min_bound`` parameters::

        sage: g = game_theory.normal_form_games.RandomImitationGame(ZZ, 5, -2, 3)
        sage: g.payoff_matrices()
        (
        [1 0 0 0 0]  [ 2 -2  0  2 -1]
        [0 1 0 0 0]  [ 2 -2 -1 -2  1]
        [0 0 1 0 0]  [-1 -2 -2  2  0]
        [0 0 0 1 0]  [-1  1  0  1  1]
        [0 0 0 0 1], [ 2  0 -1  1  2]
        )
        sage: g.obtain_nash()
        [[(0, 0, 0, 0, 1), (0, 0, 0, 0, 1)],
         [(0, 0, 0, 0, 1), (1/2, 0, 0, 0, 1/2)],
         [(0, 0, 0, 1, 0), (0, 0, 0, 1/2, 1/2)],
         [(0, 0, 0, 1, 0), (0, 0, 0, 1, 0)],
         [(0, 0, 0, 1, 0), (0, 1/3, 0, 1/3, 1/3)],
         [(0, 0, 0, 1, 0), (0, 1/2, 0, 1/2, 0)],
         [(2/11, 0, 0, 3/11, 6/11), (1/3, 0, 0, 1/3, 1/3)],
         [(1, 0, 0, 0, 0), (1/2, 0, 0, 1/2, 0)],
         [(1, 0, 0, 0, 0), (1, 0, 0, 0, 0)]]
        sage: g = game_theory.normal_form_games.RandomImitationGame(QQ, 5, 6, 8)
        sage: g.payoff_matrices()
        (
        [1 0 0 0 0]  [  -6 -3/5  1/4  3/4   -5]
        [0 1 0 0 0]  [-4/3  4/5  5/4   -3  1/3]
        [0 0 1 0 0]  [-4/7 -3/2    1 -1/4    4]
        [0 0 0 1 0]  [ 6/5  5/3  4/7  2/5  2/5]
        [0 0 0 0 1], [   1 -3/4    5    2  1/4]
        )
        sage: g.obtain_nash()
        [[(0, 0, 19/31, 0, 12/31), (0, 0, 1/2, 0, 1/2)]]

    """
    from sage.rings.all import RR, ZZ

    if ring == RR :
        if min_bound is None:
            min_bound = -1
        if max_bound is None:
            max_bound = 1

    B = matrix(ring, n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = ring.random_element(min_bound, max_bound)
    A = matrix.identity(ring, n)

    return NormalFormGame([A, B])
