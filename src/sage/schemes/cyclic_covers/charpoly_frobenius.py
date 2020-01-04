r"""

Computation of the Frobenius polynomial using Newton's identities

"""


# *****************************************************************************
#  Copyright (C) 2018 Edgar Costa <edgarc@mit.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import division
from sage.rings.integer_ring import ZZ
from sage.functions.log import log
from sage.functions.other import floor


def charpoly_frobenius(frob_matrix, charpoly_prec, p, weight, a=1, known_factor=[1]):
    """
    Return the characteristic polynomial of the given Frobenius matrix.

    INPUT:

    - ``frob_matrix`` -- a matrix representing the Frobenius matrix up to some precision

    - ``charpoly_prec`` -- a vector ai, such that, `frob_matrix.change_ring(ZZ).charpoly()[i]`
        will be correct mod `p^ai`, this can be easily deduced from the Hodge numbers and
        knowing the q-adic precision of ``frob_matrix``

    - ``p`` -- prime `p`

    - ``weight`` -- weight of the motive

    - ``a`` -- `q = q^a`

    - ``known_factor`` -- the list of coefficients of the known factor

    OUTPUT:

    A list of integers corresponding to the characteristic polynomial of the Frobenius action

    EXAMPLES::

        sage: from sage.schemes.cyclic_covers.charpoly_frobenius import charpoly_frobenius
        sage: M = Matrix([[O(17), 8 + O(17)], [O(17), 15 + O(17)]])
        sage: charpoly_frobenius(M, [2, 1, 1], 17, 1, 1)
        [17, 2, 1]

        sage: R = Zq(17**2 , names=('a',))
        sage: M = Matrix(R, [[8*17 + 16*17**2 + O(17**3), 8 + 11*17 + O(17**2)], [7*17**2 + O(17**3), 15 + 8*17 + O(17**2)]])
        sage: charpoly_frobenius(M*M, [3, 2, 2], 17, 1, 2)
        [289, 30, 1]

        sage: M = Matrix([[8*31 + 8*31**2 + O(31**3), O(31**3), O(31**3), O(31**3)], [O(31**3), 23*31 + 22*31**2 + O(31**3), O(31**3), O(31**3)], [O(31**3), O(31**3), 27 + 7*31 + O(31**3), O(31**3)], [O(31**3), O(31**3), O(31**3), 4 + 23*31 + O(31**3)]])
        sage: charpoly_frobenius(M, [4, 3, 2, 2, 2], 31, 1, 1)
        [961, 0, 46, 0, 1]

        sage: M = Matrix([(4*43^2 + O(43^3), 17*43 + 11*43^2 + O(43^3), O(43^3), O(43^3), 17 + 37*43 + O(43^3), O(43^3)),
        ....:  (30*43 + 23*43^2 + O(43^3), 5*43 + O(43^3), O(43^3), O(43^3), 3 + 38*43 + O(43^3), O(43^3)),
        ....:  (O(43^3), O(43^3), 9*43 + 32*43^2 + O(43^3), 13 + 25*43 + O(43^3), O(43^3), 17 + 18*43 + O(43^3)),
        ....:  (O(43^3), O(43^3), 22*43 + 25*43^2 + O(43^3), 11 + 24*43 + O(43^3), O(43^3), 36 + 5*43 + O(43^3)),
        ....:  (42*43 + 15*43^2 + O(43^3), 22*43 + 8*43^2 + O(43^3), O(43^3), O(43^3), 29 + 4*43 + O(43^3), O(43^3)),
        ....:  (O(43^3), O(43^3), 6*43 + 19*43^2 + O(43^3), 8 + 24*43 + O(43^3), O(43^3), 31 + 42*43 + O(43^3))])
        sage: charpoly_frobenius(M, [5, 4, 3, 2, 2, 2, 2], 43, 1, 1)
            [79507, 27735, 6579, 1258, 153, 15, 1]

        sage: M = Matrix([(1 + O(4999), O(4999), 0, 0),
        ....:  (O(4999), 4860 + O(4999), 0, 0),
        ....:  (0, 0, O(4999), O(4999)),
        ....:  (0, 0, O(4999), 1 + O(4999))])
        sage: charpoly_frobenius(M, [2, 1, 1], 4999, 1, 1, [1, -2 ,1 ])
        [4999, 139, 1]

    TESTS::

        sage: M = Matrix([[-149196156000219, 0, 0, 0, 0, 0, 0, 0],
        ....:             [0, 76324364094257, 0, 0, 0, 0, 0, 0],
        ....:             [0, 0, 76324364094257, 0, 0, 0, 0, 0],
        ....:             [0, 0, 0, -149196156000219, 0, 0, 0, 0],
        ....:             [0, 0, 0, 0, 281855171388275, 0, 0, 0],
        ....:             [0, 0, 0, 0, 0, -208983379482579, 0, 0],
        ....:             [0, 0, 0, 0, 0, 0, -208983379482579, 0],
        ....:             [0, 0, 0, 0, 0, 0, 0, 281855171388275]])
        sage: charpoly_frobenius(M, [9, 8, 7, 6, 5, 5, 5, 5, 5], 1009, 1, 2)
        [1074309286591662654798721,
         561382189105547134612,
         -2982540407204025062,
         -247015136050256,
         4390163797795,
         -242628176,
         -2877542,
         532,
         1]
        sage: M = Matrix([[0, 0, 0, -338082603, 0, 0, 0, 0],
        ....:             [0, 0, -317436968, 0, 0, 0, 0, 0],
        ....:             [0, -120741807, 0, 0, 0, 0, 0, 0],
        ....:             [200618482, 0, 0, 0, 0, 0, 0, 0],
        ....:             [0, 0, 0, 0, 0, 0, 0, 123492519],
        ....:             [0, 0, 0, 0, 0, 0, 426826171, 0],
        ....:             [0, 0, 0, 0, 0, 157417117, 0, 0],
        ....:             [0, 0, 0, 0, 373415235, 0, 0, 0]])
        sage: charpoly_frobenius(M, [7, 6, 5, 4, 3, 3, 3, 3, 3], 1009, 1, 1)
        [1036488922561, 0, 270809546, 0, -1474149, 0, 266, 0, 1]

        sage: M = Matrix({(0, 31): 1814236329200021268558465351501717,
        ....: (1, 30): 3268331092352160631300311212049390,
        ....: (2, 29): 1002349136486054751305109007707560,
        ....: (3, 28): 1789497403160078628636360424523308,
        ....: (4, 19): 919866278512654133838788268427125,
        ....: (5, 18): 2918980842679879118243999587726673,
        ....: (6, 17): 2062741569795231121341967954037400,
        ....: (7, 16): 3562554496811633214919332352788305,
        ....: (8, 7): 287823825201170974551150606916601,
        ....: (9, 6): 2657175570144838727074228404244845,
        ....: (10, 5): 3200631048273888400670606576807785,
        ....: (11, 4): 707085630754978281870563133348521,
        ....: (12, 39): 679572779843478608532167180287595,
        ....: (13, 38): 510867456922807824071915371084390,
        ....: (14, 37): 3300741705093235469798877501619286,
        ....: (15, 36): 1374430202827161695034370373469332,
        ....: (16, 27): 1897240889699239396313755822318254,
        ....: (17, 26): 3171751877741319729745976757727266,
        ....: (18, 25): 1151779650995750952707414056498421,
        ....: (19, 24): 1309748952162524211332312241346156,
        ....: (20, 15): 2914640274871541651939754878647777,
        ....: (21, 14): 2524322227034087814555116576604052,
        ....: (22, 13): 693999428630644346611319813759997,
        ....: (23, 12): 2093267437436875555592094407087011,
        ....: (24, 3): 101158112439244133585487537448909,
        ....: (25, 2): 638873050956374173808321501215560,
        ....: (26, 1): 3529335795023815426485172749287314,
        ....: (27, 0): 618726320422582798159865537548600,
        ....: (28, 35): 2510605595766272594980682702750921,
        ....: (29, 34): 2978146199632282120435531158312695,
        ....: (30, 33): 1724161588290366191539756998844438,
        ....: (31, 32): 516507426627993787229114955328811,
        ....: (32, 23): 1716672265998537901154333190869011,
        ....: (33, 22): 3787144776814278856737374038432424,
        ....: (34, 21): 3765560528316833596614887925578722,
        ....: (35, 20): 1628311006615824767735977131865996,
        ....: (36, 11): 3638935478569769465046956942756848,
        ....: (37, 10): 1878821491042105813643148323053706,
        ....: (38, 9): 1187568624951630613061547491748348,
        ....: (39, 8): 2538351040819233009959661983810741}
        ....: )
        sage: charpoly_frobenius(M,
        ....: [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
        ....:  15, 14, 13, 12] + [11]*21, 1129, 1, 1)
        [11320844849639649951608809973589776933203136765026963553258401,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         24687045654725446027864774006541463602997309796,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         20187877911930897108199045855206,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         7337188909826596,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         1]
        sage: F = matrix(Qp(19),
        ....:            [(11009078382, 709030980, 836739860, 436167534, 292239893, 788935290, 459787897, 618270108, 291407902, 813982059, 436915165, 712519076, 274108478, 842218073, 873188529, 84850808, 253606
        ....: 528, 101055490, 834667606, 392416557, 34831360),
        ....: (6703391311, 123614798, 444418569, 737529023, 823831203, 385997426, 81567418, 574826817, 76160778, 552546125, 331626722, 233773834, 124809651, 60769125, 544807026, 23912488, 765433620, 362008786,
        ....: 459200151, 832393287, 36109377),
        ....: (14027710925, 769353472, 1333496, 445807355, 565921213, 759344291, 286119575, 757107554, 254337325, 173291685, 561289431, 688306616, 802147985, 582351976, 612685115, 445508884, 627434378, 93285402
        ....: , 594225779, 237664141, 15032606),
        ....: (10230013307, 657529105, 865927622, 832543995, 25659025, 463418189, 661653017, 116594260, 476740533, 341708578, 632416615, 163609133, 729772159, 25671185, 832339536, 85964341, 649899503, 746647313
        ....: , 882364940, 475133057, 29697524),
        ....: (13544706282, 519325385, 522005031, 77374346, 593212927, 16000964, 226081627, 730667002, 864832443, 115261391, 686326208, 546089032, 279725410, 133853119, 494323608, 651859562, 309394480, 72603155
        ....: 3, 300079122, 798325185, 38356208),
        ....: (2114693958, 379215851, 2849297, 763026662, 809178574, 363345056, 516634871, 544148448, 28364967, 791976126, 312108022, 475923685, 634446290, 589784244, 526877448, 78020156, 694947534, 378288385,
        ....: 549892528, 115791776, 27178803),
        ....: (8201805924, 702420405, 740994699, 664615763, 275976558, 371434192, 598146790, 290105965, 348195672, 675097740, 132835935, 460798393, 479313988, 35743503, 391087317, 799123508, 755660343, 42443011
        ....: 3, 133110979, 504978713, 12458080),
        ....: (3785191856, 271673419, 244305002, 805279223, 680160651, 319602515, 630916470, 169264350, 122213928, 295238416, 112265908, 539354273, 790353431, 318685860, 732925760, 39166030, 127221834, 60838323
        ....: , 337119641, 870976986, 1448255),
        ....: (16749350573, 412269562, 149357480, 377515237, 843659052, 419022447, 739152212, 118583332, 82937432, 682782879, 394384824, 44425325, 103461384, 352452470, 436574476, 346452783, 182366617, 28336894
        ....: 5, 662794917, 317618782, 24381611),
        ....: (7429785403, 598980985, 357367732, 537650999, 36416996, 3551290, 492704542, 891068175, 404275920, 117359884, 261793267, 730925801, 539512391, 506466394, 741860244, 16739931, 597168613, 559371837,
        ....: 799968039, 479797709, 29444786),
        ....: (12112559717, 208269298, 464132418, 651389065, 414821186, 376370202, 42626690, 148726965, 831713068, 14123251, 611109445, 293268876, 359600004, 197934248, 422461409, 520943748, 311378517, 57729408
        ....: 1, 738958127, 515504523, 31857672),
        ....: (10118356368, 235280211, 659756247, 345436739, 256102786, 738732882, 545212543, 532594757, 338732608, 129764148, 308045746, 72149650, 558903772, 473715410, 699034852, 734871626, 326306361, 2364587
        ....: 43, 657694291, 627200317, 25892020),
        ....: (10440399775, 773503946, 828638564, 250491953, 677967063, 749991579, 523418726, 843566807, 168645748, 421663770, 819331794, 450019617, 694093902, 881515811, 255386011, 659311951, 302059587, 101610
        ....: 67, 472650365, 324078858, 33128128),
        ....: (7447735406, 61847299, 839076024, 241614165, 862220209, 863862341, 843925014, 379706488, 58299752, 8212750, 505435397, 146278188, 35720285, 710909738, 551058634, 456185516, 537754492, 545129000, 6
        ....: 44267922, 648921250, 10419599),
        ....: (13153018755, 537542129, 97194063, 395823637, 863359791, 631441497, 227091895, 771021254, 839552639, 415611567, 561558737, 594603385, 58497637, 482495633, 55060518, 460267913, 527064579, 511969915
        ....: , 346313684, 751310046, 23580385),
        ....: (3564963445, 354756562, 584228625, 301523274, 623546256, 783391508, 129133196, 58862437, 145858801, 306311749, 839640780, 242011113, 658908600, 4739911, 336057997, 10829031, 59786122, 866972812, 6
        ....: 24542293, 650575884, 42102199),
        ....: (2510113864, 40030302, 810334629, 146808649, 454276358, 236416525, 407354870, 481838195, 748232160, 696920380, 138002415, 85411726, 784262772, 408107118, 465829270, 311234288, 110620014, 150938717
        ....: , 665672087, 832237449, 36255422),
        ....: (6654231414, 882174408, 755511744, 27736751, 877076917, 578404384, 325177077, 357942672, 554557161, 453710709, 92308004, 760465291, 864441917, 12038970, 46761945, 195939666, 568849911, 584304606,
        ....: 655944486, 888687304, 24181283),
        ....: (14544424304, 826543206, 796282628, 23338517, 824298603, 302094148, 660109172, 718095101, 706279229, 526516847, 141230895, 790244903, 856267452, 573246036, 211086922, 372229190, 691861155, 4199917
        ....: 51, 436860141, 622045522, 8325584),
        ....: (15767321551, 661362716, 507371364, 195363985, 443294472, 225547499, 82193848, 636150685, 403093588, 90068531, 12901285, 353476874, 339244848, 517694292, 69606861, 789231082, 26751924, 494426075,
        ....: 103926086, 745554015, 41310692),
        ....: (16927354258, 644384221, 804662711, 109546134, 859428330, 531587396, 742764321, 885729859, 218211295, 554897109, 654169354, 771389569, 660454516, 64500193, 685951224, 96204372, 320976804, 52781829
        ....: , 293886376, 814843785, 2152625)]
        ....: )
        ....: F+= F.base_ring()(0).add_bigoh(6)*ones_matrix(*F.dimensions())
        ....: charpoly_frobenius(F, [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 10, 8, 7, 6, 6, 6], 19
        ....: , 2)
        [-714209495693373205673756419,
         265108233858482020942613186,
         -42900468322739431511530492,
         4082888209049677628524555,
         -270269604512861034117097,
         14391708424321308787452,
         -647994422169919022131,
         21573180516137871267,
         -294370884235799413,
         -21783678413966953,
         2262363449128569,
         -119071760480451,
         3175926288667,
         118884941287,
         -24134536953,
         2008116289,
         -123544308,
         6426883,
         -268945,
         7828,
         -134,
         1]


    """
    assert known_factor[-1] == 1
    try:
        cp = frob_matrix.change_ring(ZZ).charpoly().list()
    except ValueError:
        # the given matrix wasn't integral
        cp = frob_matrix.charpoly().change_ring(ZZ).list()
    assert len(charpoly_prec) == len(cp) - (len(known_factor) - 1)
    assert cp[-1] == 1

    # reduce cp mod prec
    degree = len(charpoly_prec) - 1
    halfdegree = floor(degree / 2) + 1
    mod = [0] * (degree + 1)
    for i in range(len(charpoly_prec)):
        mod[-i] = p**charpoly_prec[-i]
        cp[-i] = cp[-i] % mod[-i]

    # figure out the sign
    if weight % 2 == 1:
        # for odd weight the sign is always 1
        # it's the charpoly of a USp matrix
        # and charpoly of a symplectic matrix is reciprocal
        sign = 1
    else:
        # For the moment I will not worry about this case
        if known_factor != [1]:
            raise NotImplementedError()
        for i in range(halfdegree):
            p_power = p**min(
                charpoly_prec[i],
                charpoly_prec[degree - i] + (a * (degree - 2 * i) * weight / 2),
            )
            # Note: degree*weight = 0 mod 2
            if cp[i] % p_power != 0 and cp[degree - i] % p_power != 0:
                other = cp[degree - i] * p**(a * (degree - 2 * i) * weight / 2)
                if (cp[i] + other) % p_power == 0:
                    sign = -1
                else:
                    sign = 1
                assert (-sign * cp[i] + other) % p_power == 0
                break
    cp[0] = sign * p**(a * degree * weight / 2)

    # calculate the i-th power sum of the roots and correct cp along the way
    e = cp[-halfdegree:]
    e.reverse()
    for k in range(halfdegree):
        if k % 2 != 0:
            e[k] = -e[k] % mod[degree - k]
        # e[k] = cp[degree - k] if (k%2 ==0) else -cp[degree - k]
        if k > 0:
            # verify if p^charpoly_prec[degree - k] > 2*degree/k * q^(w*k/2)
            assert (
                log(k) / log(p) + charpoly_prec[degree - k]
                > log(2 * degree) / log(p) + a * 0.5 * weight * k
            ), (
                "log(k)/log(p) + charpoly_prec[degree - k] <= log(2*degree)/log(p) + a*0.5*weight*k, k = %d"
                % k
            )

    fix_e = known_factor[:]
    fix_e.reverse()
    if len(fix_e) < halfdegree:
        fix_e.extend([0] * (halfdegree - len(fix_e)))
    for i in range(halfdegree):
        if i % 2 != 0:
            fix_e[i] *= -1

    # e[k] = \sum x_{i_1} x_{i_2} ... x_{i_k} # where x_* are eigenvalues
    # and i_1 < i_2 ... < i_k

    # s[k] = \sum x_i ^k for k>0
    s = [None] * (halfdegree)
    res = [None] * len(charpoly_prec)
    res[0] = sign * p**(a * degree * weight / 2)
    res[-1] = 1
    e[1] -= fix_e[1]
    e[1] = e[1] % mod[degree - 1]
    for k in range(1, halfdegree):
        # assume that s[i] and e[i] are correct for i < k
        # e[k] correct modulo mod[degree - k]
        # S = sum (-1)^i e[k-i] * s[i]
        # s[k] = (-1)^(k-1) (k*e[k] + S) ==> (-1)^(k-1) s[k] - S = k*e[k]
        S = sum((-1)**i * e[k - i] * s[i] for i in range(1, k))
        s[k] = (-1)**(k - 1) * (S + k * e[k])
        # hence s[k] is correct modulo k*mod[degree - k]
        localmod = k * mod[degree - k]
        # s[k] +=   (-1)**k * fix_power_sum[k]
        s[k] = s[k] % localmod

        # |x_i| = p^(w*0.5)
        # => s[k] <= degree*p^(a*w*k*0.5)
        # recall, 2*degree*p^(a*w*k*0.5) /k < mod[degree - k]
        if s[k]**2 > degree**2 * p**(a * weight * k):
            s[k] = -(-s[k] % localmod)

        # now correct e[k] with:
        # (-1)^(k-1) s[k] - S = k*e[k]
        e[k] = (-S + (-1)**(k - 1) * s[k]) // k
        assert (-S + (-1)**(k - 1) * s[k]) % k == 0
        res[degree - k] = e[k] if k % 2 == 0 else -e[k]
        res[k] = sign * res[degree - k] * p**(a * (degree - 2 * k) * weight / 2)
        # fix e[k + 1]
        if k + 1 < halfdegree:
            e[k + 1] -= sum([fix_e[k + 1 - i] * e[i] for i in range(k + 1)])
            e[k + 1] = e[k + 1] % mod[degree - (k + 1)]
    return res
