r"""
Poisson tensors

AUTHORS:

- Tobias Diez (2020): initial version
"""

# *****************************************************************************
#       Copyright (C) 2020 Tobias Diez
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************


from sage.manifolds.differentiable.multivectorfield import (
    MultivectorField,
    MultivectorFieldParal,
)
from sage.manifolds.differentiable.diff_form import DiffForm
from sage.manifolds.differentiable.vectorfield import VectorField
from sage.manifolds.differentiable.scalarfield import DiffScalarField
from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from typing import Optional, Union


class PoissonTensorField(MultivectorField):
    r"""
    A Poisson bivector field `\varpi` on a differentiable manifold.

    That is, at each point `m \in M`, `\varpi_m` is a bilinear map of the type:

    .. MATH::

        \varpi_m:\ T^*_m M \times T^*_m M  \to \RR

    where `T^*_m M` stands for the cotangent space to the
    manifold `M` at the point `m`, such that `\varpi_m` is skew-symmetric and the
    Schouten bracket of `\varpi` with itself vanishes.
    """

    def __init__(
        self,
        manifold: Union[DifferentiableManifold, VectorFieldModule],
        name: Optional[str] = "varpi",
        latex_name: Optional[str] = "\\varpi",
    ):
        r"""
        Construct a Poisson bivector field.

        INPUT:

        - ``manifold`` -- module `\mathfrak{X}(M)` of vector
        fields on the manifold `M`, or the manifold `M` itself
        - ``name`` -- (default: ``varpi``) name given to the Poisson tensor
        - ``latex_name`` -- (default: ``\\varpi``) LaTeX symbol to denote the Poisson tensor;
        if ``None``, it is formed from ``name``

        EXAMPLES:

        Standard Poisson tensor on `\RR^2`::

            sage: from sage.manifolds.differentiable.poisson_tensor import PoissonTensorField
            sage: M.<q, p> = EuclideanSpace(2, "R2", symbols=r"q:q p:p")
            sage: poisson = PoissonTensorField(M, 'varpi'); poisson
            2-vector field varpi on the Euclidean plane R2
        """
        try:
            vector_field_module = manifold.vector_field_module()
        except AttributeError:
            vector_field_module = manifold

        MultivectorField.__init__(
            self, vector_field_module, 2, name=name, latex_name=latex_name
        )

    def hamiltonian_vector_field(self, function: DiffScalarField) -> VectorField:
        r"""
        Return the Hamiltonian vector field `X_f` generated by the given function `f: M \to \RR`.

        The Hamiltonian vector field is defined by
        .. MATH::

            X_f = - \varpi^\sharp (d f),

        where `\varpi^\sharp: T^* M \to TM` is given by
        `\beta(\varpi^\sharp(\alpha)) = \varpi(\alpha, \beta)`.

        INPUT:

        - ``function`` -- the function generating the Hamiltonian vector field

        EXAMPLES:

            sage: M.<q, p> = EuclideanSpace(2, "R2", symbols=r"q:q p:p")
            sage: poisson = M.poisson_tensor('varpi')
            sage: poisson.set_comp()[1,2] = -1
            sage: f = M.scalar_field({ chart: function('f')(*chart[:]) for chart in M.atlas() }, name='f')
            sage: Xf = poisson.hamiltonian_vector_field(f)
            sage: Xf.display()
            Xf = d(f)/dp e_q - d(f)/dq e_p
        """
        vector_field = -self.sharp(function.exterior_derivative())
        if function._name is not None:
            vector_field.set_name(f"X{function._name}", f"X_{{{function._latex_name}}}")
        return vector_field

    def sharp(self, form: DiffForm) -> VectorField:
        r"""
        Return the image of the given differential form under the map `\varpi^\sharp: T^* M \to TM` defined by
        .. MATH::
            `\beta(\varpi^\sharp(\alpha)) = \varpi(\alpha, \beta)`.

        for all `\alpha, \beta \in T^*_m M`.

        In indicies, `\alpha^i = \varpi^{ij} \alpha_j`.

        INPUT:

        - ``form`` -- the differential form to calculate it's sharp of

        EXAMPLES:

            sage: M.<q, p> = EuclideanSpace(2, "R2", symbols=r"q:q p:p")
            sage: poisson = M.poisson_tensor('varpi')
            sage: poisson.set_comp()[1,2] = -1
            sage: a = M.one_form(1, 0, name='a')
            sage: poisson.sharp(a).display()
            a_sharp = e_p
        """
        if form.degree() != 1:
            raise ValueError(
                f"the degree of the differential form must be one but it is {form.degree()}"
            )

        vector_field = form.up(self)
        vector_field.set_name(f"{form._name}_sharp", f"{form._latex_name}^\\sharp")
        return vector_field

    def poisson_bracket(
        self, f: DiffScalarField, g: DiffScalarField
    ) -> DiffScalarField:
        r"""
        Return the Poissen bracket
        .. MATH::
            {f, g} = \varpi(\dif f, \dif g)

        of the given functions.

        INPUT:

        - ``f`` -- first function
        - ``g`` -- second function

        EXAMPLES:

            sage: M.<q, p> = EuclideanSpace(2, "R2", symbols=r"q:q p:p")
            sage: poisson = M.poisson_tensor('varpi')
            sage: poisson.set_comp()[1,2] = -1
            sage: f = M.scalar_field({ chart: function('f')(*chart[:]) for chart in M.atlas() }, name='f')
            sage: g = M.scalar_field({ chart: function('g')(*chart[:]) for chart in M.atlas() }, name='g')
            sage: poisson.poisson_bracket(f, g).display()
            poisson(f, g): R2 → ℝ
               (q, p) ↦ d(f)/dp*d(g)/dq - d(f)/dq*d(g)/dp
        """
        poisson_bracket = self.contract(0, f.exterior_derivative()).contract(
            0, g.exterior_derivative()
        )
        poisson_bracket.set_name(
            f"poisson({f._name}, {g._name})",
            "\\{" + f"{f._latex_name}, {g._latex_name}" + "\\}",
        )
        return poisson_bracket


class PoissonTensorFieldParal(PoissonTensorField, MultivectorFieldParal):
    def __init__(
        self,
        manifold: Union[DifferentiableManifold, VectorFieldModule],
        name: Optional[str] = None,
        latex_name: Optional[str] = None,
    ):
        PoissonTensorField.__init__(self, manifold, name, latex_name)
        MultivectorFieldParal.__init__(self, self._vmodule, 2, name, latex_name)
