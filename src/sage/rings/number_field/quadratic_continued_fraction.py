##########################################
#
# Authors: Xevi Guitart and Marc Masdeu
#
##########################################

from collections import deque, namedtuple
from sage.functions.all import floor, ceil
from sage.matrix.constructor import Matrix
from sage.misc.all import union
from sage.plot import plot, polygon, point
from sage.rings.all import ZZ, QQ, RR
from sage.structure.sage_object import SageObject


class Regions(SageObject):
    def __init__(self, parent):
        self._regions = dict([])
        self._F = parent._F
        self._eps = parent._eps
        self._epsto = [self._F(1)]
        self.fundom_rep = parent.fundom_rep
        self.embed = parent.embed
        self._parent = parent
        self._N = 1  # Means that we have regions up to _N-1 = 0, so none
        self._area = self._parent._Rw

    def get_regions(self, B):
        regions = self._regions
        if B > self._N:
            eps = self._eps
            newelts = {(ZZ(nv[0]), [I.gens_reduced()[0]
                                    for I in nv[1] if I.is_principal()])
                       for nv in self._F.ideals_of_bdd_norm(B - 1).iteritems()
                       if nv[0] >= self._N}
            self._N = B
            for nn in range(len(self._epsto), B):
                self._epsto.append(self._epsto[nn - 1] * eps)
            for nn, v in newelts.iteritems():
                assert not(nn in regions)
                regions[nn] = []
                for q in v:
                    invden = 1 / q
                    invden_red = self.fundom_rep(invden)
                    q0 = invden_red - invden
                    m_invden_red = self.fundom_rep(-invden)
                    m_q0 = m_invden_red + invden

                    a, b = self._parent._change_basis(invden_red)
                    am, bm = self._parent._change_basis(m_invden_red)

                    regions[nn].extend([[invden_red, [q0, q],
                                         Region(self.embed(invden_red),
                                                RR(1 / nn)),
                                         [a, b], nn],
                                        [m_invden_red, [m_q0, -q],
                                         Region(self.embed(m_invden_red),
                                                RR(1 / nn)),
                                         [am, bm], nn]])

                    for jj in range(1, nn):
                        x = self._epsto[jj] * invden
                        x_red = self.fundom_rep(x)
                        x_red_minus = self.fundom_rep(-x_red)
                        regs = regions[nn]
                        if x_red == regs[0][0] or x_red_minus == regs[0][0]:
                            continue

                        q0_minus = x_red_minus + x
                        q0 = x_red - x
                        q1 = 1 / x
                        a, b = self._parent._change_basis(x_red)
                        am, bm = self._parent._change_basis(x_red_minus)
                        regions[nn].extend([[x_red, [q0, q1],
                                             Region(self.embed(x_red),
                                                    RR(1 / nn)),
                                             [a, b], nn],
                                            [x_red_minus, [q0_minus, -q1],
                                             Region(self.embed(x_red_minus),
                                                    RR(1 / nn)),
                                             [am, bm], nn]])

        else:
            # We have all the required regions
            pass
        return [reg for nn, lst in regions.iteritems()
                if nn < B for reg in lst]


class Hole(SageObject):
    def __init__(self, parent, v, depth=0):
        self._parent = parent
        if isinstance(v, list):
            #Input in as [xmin,ymin,xmax,ymax]
            if not(v[0] < v[2] and v[1] < v[3]):
                raise ValueError('The input should be of the form '
                                 '[xmin, ymin, xmax, ymax], with '
                                 'xmin < xmax and ymin < ymax.')
            (self.xmin, self.ymin, self.xmax, self.ymax) = v
        elif isinstance(v, Hole):
            self.xmin = v.xmin
            self.ymin = v.ymin
            self.xmax = v.xmax
            self.ymax = v.ymax
        self.size = max([self.xmax - self.xmin, self.ymax - self.ymin])
        self._corners = None
        self._depth = depth

    def depth(self):
        return self._depth

    def corners(self):
        if self._corners is None:
            Pt = self._parent.Pointxy
            self._corners = [Pt(self.xmin, self.ymin),
                             Pt(self.xmin, self.ymax),
                             Pt(self.xmax, self.ymax),
                             Pt(self.xmax, self.ymin)]
        return self._corners

    def plot(self, color='blue'):
        return plot(polygon(self.corners()))

    def overlaps_fundamental_domain(self):
        return any([self._parent.in_fundamental_domain(x)
                    for x in self.corners()])

    def _repr_(self):
        return 'Hole with corners = ({},{}),({},{})'.format(self.xmin,
                                                            self.ymin,
                                                            self.xmax,
                                                            self.ymax)


class Region(SageObject):
    def __init__(self, center, radius):
        #assert(isinstance(center, Pointxy))
        self._center = center
        radius = abs(radius)
        if radius > 1:
            print radius
            print center
            assert(0)
        self._radius = radius

    def radius(self):
        return self._radius

    def center(self):
        return self._center

    def plot_center(self, color='red'):
        return plot(point(self._center, color=color))

    def plot(self, color='red'):
        r = 2 * self._radius
        a = self._center[0]
        b = self._center[1]
        Pts = [(a, b + r), (a + r, b), (a, b - r), (a - r, b)]
        return plot(polygon(Pts), color=color)

    def contains_point(self, P):
        return (abs((P.x - self._center.x) * (P.y - self._center.y)) <
                self._radius)

    def contains_hole(self, h):
        return all([self.contains_point(P) for P in h.corners()])


class QuadraticContinuedFraction(SageObject):
    def __init__(self, F, Nbound=50, Tbound=5):
        self._F = F
        self._solved = False
        self._disc = self._F.discriminant()
        self._w = self._F.maximal_order().ring_generators()[0]
        self._eps = self._F.units()[0]
        self._Phi = self._F.real_embeddings(prec=300)

        self.Pointxy = namedtuple('Pointxy', 'x y')
        if self._disc % 2 == 0:
            self._changebasismatrix = 1
        else:
            self._changebasismatrix = Matrix(ZZ, 2, 2, [1, -1, 0, 2])

        # We assume here that Phi orders the embeddings to that
        # Phi[1] is the largest
        self._Rw = self._Phi[1](self._w)
        self._Rwconj = self._Phi[0](self._w)

        self._used_regions = []
        self._Nboundmin = 2
        self._Nboundmax = Nbound
        self._Nboundinc = 1
        self._Nbound = Nbound
        self._epsto = [self._F(1)]
        self._Dmax = self.embed(self._w)
        self._Tbound = Tbound
        self._ranget = sorted(range(-self._Tbound, self._Tbound + 2), key=abs)
        self._M = Matrix(RR, 2, 2, [1, self._Rw, 1, self._Rwconj]).inverse()
        self.initialize_fundom()
        self._master_regs = Regions(self)
        self._maxdepth = 0

    def fundom_rep(self, x):
        """
        Takes an element in the number field F and returns
        y in the fundamental domain such that such that (x - y) belongs to R
        """
        # Now xp contains a,b such that x=a+b sqrt(D)
        y00, y10 = self._change_basis(x)
        d = y00.floor() + self._w * y10.floor()
        return x - d

    def _change_basis(self, x):
        v = x.parts()
        tmp = self._changebasismatrix * Matrix(QQ, 2, 1, v)
        return tmp[0, 0], tmp[1, 0]

    def rangea_gen(self, xmin, ymin, xmax, ymax):
        """
        Return a function that, given a t-value, will return a range of
        a-values.

        The arguments specify the region that we want to cover.
        """
        wprime0 = self._Rw
        wprime1 = self._Rwconj

        def rangea(t):
            tp0 = t * wprime0
            tp1 = t * wprime1
            return sorted(list(union(range(floor(xmin - 1 - tp0 - wprime0),
                                           ceil(xmax - tp0) + 1),
                                     range(floor(ymin - 1 - tp1 - wprime1),
                                           ceil(ymax - tp1) + 1))), key=abs)
        return rangea

    def in_fundamental_domain(self, P):
        #assert(isinstance(P,Pointxy))
        A = self._M * Matrix(RR, 2, 1, [P.x, P.y])
        return A[0, 0] >= 0 and A[0, 0] < 1 and A[1, 0] >= 0 and A[1, 0] < 1

    def embed(self, x):
        """
        Given x in F (or in R) embeds it into R^2
        """
        return self.Pointxy(self._Phi[1](x), self._Phi[0](x))

    def embed_coords(self, a, b):
        return self.Pointxy(a + b * self._Rw, a + b * self._Rwconj)
        #return ((ZZ(a)+self._Rw*ZZ(b),ZZ(a)+self._Rwconj*ZZ(b)))

    def initialize_fundom(self):
        """
        Build the fundamental domain
        """
        Pts = [self.embed_coords(a, b) for a in [0, 1] for b in [0, 1]]
        xmin = min([P.x for P in Pts])
        xmax = max([P.x for P in Pts])
        ymin = min([P.y for P in Pts])
        ymax = max([P.y for P in Pts])

        self._aplusbomega = dict([])
        rangea = self.rangea_gen(xmin - 1, ymin - 1, xmax + 1, ymax + 1)
        for b in self._ranget:
            bw = b * self._w
            for a in rangea(b):
                self._aplusbomega[(a, b)] = a + bw

        dx = RR(0.49)
        dy = dx
        Nx = ceil((xmax - xmin) / dx)
        Ny = ceil((ymax - ymin) / dy)
        self._holes = deque([])
        for ii in range(Nx):
            for jj in range(Ny):
                self._holes.append(Hole(self, [xmin + ii * dx,
                                               ymin + jj * dy,
                                               xmin + (ii + 1) * dx,
                                               ymin + (jj + 1) * dy]))

    def plot_holes(self):
        myplot = plot([])
        for h in self._holes:
            myplot += h.plot()
        return myplot

    def remaining_holes(self):
        return len(self._holes)

    def test_regions(self, corners, at_least_one, embed_coords,
                     regions, Region, aplusbomega, ranget, rangea):
        ii = 0
        NumRegions = len(regions)
        #print "Looking at ",NumRegions," regions..."
        atws = [[aplusbomega[(a, b)], embed_coords(a, b), [a, b]]
                for b in ranget for a in rangea(b)]
        for ii in range(NumRegions):
            reg = regions[ii]
            for inc in atws:
                passed = [False] * 4
                a = inc[1].x
                b = inc[1].y
                for ii, P in enumerate(corners):
                    if reg[2].contains_point(self.Pointxy(P.x - a, P.y - b)):
                        at_least_one[ii] = True
                        passed[ii] = True
                    else:
                        break
                if all(passed):
                    rparts = reg[3]
                    incparts = inc[2]
                    return [reg[0] + inc[0],
                            [reg[1][0] + inc[0], reg[1][1]],
                            Region(embed_coords(rparts[0] + incparts[0],
                                                rparts[1] + incparts[1]),
                                   reg[2].radius()), reg[3], reg[4]]
        return None

    def evaluate_hole(self, h, verifying=False, maxdepth=-1):
        """
        # Given a hole, returns:
        #  -2 if at least one of the points is not covered by the regions
        #  -1 if no region covers all the points
        #  0 if the a region can cover all the points
        """
        if not h.overlaps_fundamental_domain():
            return 0
        used = self._used_regions
        corners = h.corners()
        if verifying:
            at_least_one = [False] * 4

        for reg in used:
            passed = True
            for ii, P in enumerate(corners):
                if not reg[2].contains_point(P):
                    passed = False
                    if not verifying:
                        break
                elif verifying:
                    at_least_one[ii] = True
            if passed:
                return 0

        if verifying:
            if all(at_least_one):
                return -1
            else:
                if h.depth() > maxdepth:
                    return -2
                else:
                    return -1

        rangea = self.rangea_gen(h.xmin, h.ymin, h.xmax, h.ymax)
        at_least_one = [False] * 4

        l = 100
        B = ceil(4 / (l * (h.size ** 2)))
        B = max([self._Nboundmin, B])
        #print 'B=',B,'Size ~ 10^',RealField(prec=10)((log(h.size)/log(10)))
        regions = self._master_regs.get_regions(min([B, self._Nbound]))
        nreg = self.test_regions(corners, at_least_one,
                                 embed_coords=self.embed_coords,
                                 regions=regions, Region=Region,
                                 aplusbomega=self._aplusbomega,
                                 ranget=self._ranget, rangea=rangea)
        if not(nreg is None):
            used.append(nreg)
            return 0

        return -1

    def evaluate_pointxy(self, x, y):
        P = self.Pointxy(x, y)
        for reg in self._used_regions:
            if reg[2].contains_point(P):
                return reg
        return -1

    def evaluate_number(self, x):
        y = self.fundom_rep(x)
        d = x - y
        P = self.embed(y)
        for reg in self._used_regions:
            if(reg[2].contains_point(P)):
                return [reg[1][0] + d, reg[1][1]]
        return -1

    def verify(self):
        maxdepth = self._maxdepth
        self.initialize_fundom()
        percent = RealField(prec=20)(0)
        holes = self._holes
        self._minholes = len(holes)
        while holes:
            if len(holes) - self._minholes > 1000:
                assert(0)
            h = holes.popleft()
            result = self.evaluate_hole(h, verifying=True, maxdepth=maxdepth)
            if result == -2:
                print "The solution is incorrect."
                print "The {} is not covered by any region.".format(h)
                return False
            if result == -1:
                xh = (h.xmax + h.xmin) / 2
                yh = (h.ymax + h.ymin) / 2
                depth = h.depth() + 1
                self._maxdepth = max([self._maxdepth, depth])
                h1 = Hole(self, [h.xmin, h.ymin, xh, yh], depth=depth)
                h2 = Hole(self, [h.xmin, yh, xh, h.ymax], depth=depth)
                h3 = Hole(self, [xh, h.ymin, h.xmax, yh], depth=depth)
                h4 = Hole(self, [xh, yh, h.xmax, h.ymax], depth=depth)
                holes.extendleft([h1, h2, h3, h4])
            else:
                percent += 100 * (4 ** (-h.depth()))
                #print 'Finished hole ',self._minholes,' up to ',percent,'%'

            if len(holes) < self._minholes:
                self._minholes = len(holes)
                percent = RealField(prec=20)(0)
                #print "Remaining ",self._minholes,' holes.'
        return True

    def solve(self, verify=False):
        percent = RealField(prec=20)(0)
        holes = self._holes
        self._minholes = len(holes)
        while holes:
            if len(holes) - self._minholes > 1000:
                assert(0)
            h = holes.popleft()
            result = self.evaluate_hole(h)
            if result == -2:
                self._Nboundmax = max([self._Nboundmax, self._Nbound])
                self._holes.append(h)
            if result == -1:
                xh = (h.xmax + h.xmin) / 2
                yh = (h.ymax + h.ymin) / 2
                depth = h.depth() + 1
                self._maxdepth = max([self._maxdepth, depth])
                h1 = Hole(self, [h.xmin, h.ymin, xh, yh], depth=depth)
                h2 = Hole(self, [h.xmin, yh, xh, h.ymax], depth=depth)
                h3 = Hole(self, [xh, h.ymin, h.xmax, yh], depth=depth)
                h4 = Hole(self, [xh, yh, h.xmax, h.ymax], depth=depth)

                holes.extendleft([h1, h2, h3, h4])
            else:
                percent += 100 * (4 ** (-h.depth()))
                #print 'Finished hole ',self._minholes,' up to ',percent,'%'

            if len(holes) < self._minholes:
                self._minholes = len(holes)
                percent = RealField(prec=20)(0)
                #print "Remaining ",self._minholes,' holes.'
        if verify:
            assert(self.verify())
        self._solved = True
        return [self._disc, self._maxdepth,
                1 + max([abs(self._F(reg[1][1]).norm())
                         for reg in self._used_regions]),
                len(self._used_regions), self._used_regions]
