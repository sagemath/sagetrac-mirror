import sage.groups.group as group


class GSp4_Arithmetic_Subgroups(group.Group):
    """
    Arithmetic subgroups (finite index subgroups of {GSp}(4,QQ)).

    EXAMPLES::

        sage: sage.modular.siegel.siegel_modular_group.GSp4_Arithmetic_Subgroups()
        Generic arithmetic subgroup of GSp(4,Q)
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: sage.modular.siegel.siegel_modular_group.GSp4_Arithmetic_Subgroups()
            Generic arithmetic subgroup of GSp(4,Q)
        """
        return "Generic arithmetic subgroup of GSp(4,Q)"

    def __reduce__(self):
        """
        Helper for pickle
        """
        raise NotImplementedError("all subclasses must define a __reduce__ method")


class Sp4Z_class(GSp4_Arithmetic_Subgroups):
    """
    The full Siegel modular group {Sp}(4,ZZ), regarded as a congruence subgroup of itself.

    EXAMPLES::

        sage: G = Sp4Z; G
        Siegel Modular Group Sp(4,Z)

        sage: G.gens()
        Traceback (most recent call last):
        ...
        NotImplementedError: Number of generators not known.
    """
    def level(self):
        return 1

    def _repr_(self):
        """
        Return the string representation of self
        """
        return "Siegel Modular Group Sp(4,Z)"

    def __reduce__(self):
        return _Sp4Z_ref, ()


def _Sp4Z_ref():
    """
    Return Sp4Z, (Used for pickling Sp4Z)

    EXAMPLES::

        sage: sage.modular.siegel.siegel_modular_group._Sp4Z_ref() is Sp4Z
        True
    """
    return Sp4Z


Sp4Z = Sp4Z_class()

_sp4z_gamma0_cache = {}


def Sp4Z_Gamma0_constructor(N):
    """
    Return the congruence subgroup Gamma0(N) of the Siegel modular group.

    EXAMPLES::

        sage: G = Sp4Z_Gamma0(51) ; G # indirect doctest
        Siegel Modular Group Gamma0(51)
        sage: G == Sp4Z_Gamma0(51)
        True
    """
    if N == 1:
        return Sp4Z
    try:
        return _sp4z_gamma0_cache[N]
    except KeyError:
        _sp4z_gamma0_cache[N] = Sp4Z_Gamma0_class(N)
        return _sp4z_gamma0_cache[N]


class Sp4Z_Gamma0_class(GSp4_Arithmetic_Subgroups):
    """
    The congruence subgroup Gamma_0(N) of the Siegel modular group

    EXAMPLES::

        sage: G = Sp4Z_Gamma0(11); G
        Siegel Modular Group Gamma0(11)
        sage: loads(G.dumps()) == G
        True
    """

    def __init__(self, level):
        self.__level = level

    def __repr__(self):
        return "Siegel Modular Group Gamma0(%d)" % self.level()

    def level(self):
        return self.__level

    def __reduce__(self):
        return Sp4Z_Gamma0_constructor, (self.level(),)
