TESTS::

    P.<x,y> = ProjectiveSpace(QQ, 1)
    F = DynamicalSystem_projective([x^4, y^4])
    F.is_chebyshev()


    P.<x,y> = ProjectiveSpace(QQ, 1)
    F = DynamicalSystem_projective([x^2 + y^2, y^2])
    F.is_chebyshev()


    P.<x,y> = ProjectiveSpace(QQ, 1)
    F = DynamicalSystem_projective([2*x^2 - y^2, y^2])
    F.is_chebyshev()


    P.<x,y> = ProjectiveSpace(QQ, 1)
    F = DynamicalSystem_projective([x^3, 4*y^3 - 3*x^2*y])
    F.is_chebyshev()

