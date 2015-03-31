load('padicperiods.sage')
set_verbose(1)
# Set a particular seed
# magma.eval("SetSeed(1251311619)")
# mseed = ZZ.random_element(10000)
mseed = 8033
magma.eval("SetSeed(%s)"%mseed)
print 'mseed = %s'%mseed

# A good curve
# Conductor 357 = 3 * 7 * 17
x = QQ['x'].gen()
Crv = HyperellipticCurve(x^6+8*x^4-8*x^3+20*x^2-12*x+12)
I2g, I4g, I6g, I10g = Crv.igusa_clebsch_invariants()
j1g = I2g**5 / I10g # One of absolute's Igusa invariants
j2g = I2g**3 * I4g / I10g
j3g = I2g**2 * I6g / I10g # One of absolute's Igusa invariants

p = 7
D = 3*17
sign = 1
prec = 70
working_prec = 200
x = QQ['x'].gen()
pol = x^2 - 2

path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *
G = BigArithGroup(p,D,1)
Coh = CohomologyGroup(G.Gpn)
q = ZZ(2)
K0 = (Coh.involution_at_infinity_matrix()-sign).right_kernel().change_ring(QQ)
Aq0 = Coh.hecke_matrix(q,use_magma = True).transpose().change_ring(QQ)
good_component = None
for U0,is_irred in Aq0.decomposition_of_subspace(K0):
    if U0.dimension() == 2 and is_irred and Aq0.restrict(U0).charpoly() == pol:
        print Aq0.restrict(U0).charpoly()
        assert good_component is None
        good_component = U0
assert good_component is not None
print good_component

flist = []
for row0 in (good_component.denominator()*good_component.matrix()).rows():
    col0 = [ZZ(o) for o in row0.list()]
    f = sum([a*Coh.gen(i) for i,a in enumerate(col0) if a != 0],Coh(0))
    flist.append(f)

wp = G.wp()
B = G.Gpn.abelianization()
C = G.Gn.abelianization()
Bab = B.abelian_group()
Cab = C.abelian_group()
verbose('Finding f...')
fdata = [B.ab_to_G(o).quaternion_rep for o in B.gens_small()]
# verbose('fdata = %s'%fdata)
f = B.hom_from_image_of_gens_small([C.G_to_ab(G.Gn(o)) for o in fdata])
verbose('Finding g...')
gdata = [wp**-1 * o * wp for o in fdata]
# verbose('gdata = %s'%gdata)
g = B.hom_from_image_of_gens_small([C.G_to_ab(G.Gn(o)) for o in gdata])
fg = direct_sum_of_maps([f,g])
V = Bab.gen(0).lift().parent()
good_ker = V.span_of_basis([o.lift() for o in fg.kernel().gens()]).LLL().rows()
ker = [B.ab_to_G(Bab(o)).quaternion_rep for o in good_ker]

foundg0 = False
foundg1 = False
for o in ker:
    a,b = flist[0].evaluate(o)[0], flist[1].evaluate(o)[0]
    if not foundg0 and b == 0 and a != 0:
        foundg0 = True
        g0 = o
        g0coeff = a
    if not foundg1 and a == 0 and b != 0:
        foundg1 = True
        g1 = o
        g1coeff = b
assert foundg0 and foundg1
c = lcm(g0coeff,g1coeff)
g0 = g0**ZZ(c/g0coeff)
g1 = g1**ZZ(c/g1coeff)

print [flist[0].evaluate(g0),flist[0].evaluate(g1)]
print [flist[1].evaluate(g0),flist[1].evaluate(g1)]

from homology import *
xi10,xi20 = lattice_homology_cycle(G,G.Gn(g0),G.Gn(wp**-1 * g0 * wp),working_prec)
xi11,xi21 = lattice_homology_cycle(G,G.Gn(g1),G.Gn(wp**-1 * g1 * wp),working_prec)
Phif = get_overconvergent_class_quaternionic(p,flist[0],G,prec,sign,progress_bar = True)
Phig = get_overconvergent_class_quaternionic(p,flist[1],G,prec,sign,progress_bar = True)

from integrals import *
num = integrate_H1(G,xi10,Phif,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
den = integrate_H1(G,xi20,Phif,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
A = num/den

num = integrate_H1(G,xi11,Phif,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
den = integrate_H1(G,xi21,Phif,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
B = num/den

num = integrate_H1(G,xi11,Phig,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
den = integrate_H1(G,xi21,Phig,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
D = num/den

A = (A + O(p**(prec+A.valuation()))).trace()/A.parent().degree()
B = (B + O(p**(prec+B.valuation()))).trace()/B.parent().degree()
D = (D + O(p**(prec+D.valuation()))).trace()/D.parent().degree()

F = A.parent()
T = Aq0.restrict(good_component)
TF = T.change_ring(F)
a,b = p_adic_l_invariant(A,B,D,TF)
a = a.trace()/a.parent().degree()
b = b.trace()/b.parent().degree()

# Below are precomputed values
a = 4*7 + 7^2 + 5*7^3 + 4*7^4 + 7^5 + 4*7^7 + 2*7^8 + 7^11 + 4*7^13 + 2*7^14 + 5*7^15 + 6*7^16 + 2*7^17 + 4*7^18 + 4*7^19 + 5*7^20 + 7^21 + 4*7^22 + 6*7^23 + 2*7^25 + 6*7^26 + 6*7^27 + 2*7^28 + 7^31 + 4*7^32 + 3*7^33 + 3*7^34 + 7^35 + 6*7^37 + 7^38 + 3*7^39 + 5*7^40 + 7^41 + 2*7^42 + 7^43 + 3*7^44 + 7^45 + 7^46 + 3*7^47 + 7^49 + 6*7^50 + 5*7^51 + 5*7^52 + 6*7^53 + 5*7^55 + 3*7^56 + 5*7^57 + 4*7^58 + 3*7^59 + 3*7^60 + 4*7^61 + 3*7^62 + 7^63 + 3*7^64 + 4*7^65 + 2*7^67 + 5*7^68 + 4*7^69 + 6*7^70 + 5*7^71 + 6*7^73 + 6*7^74 + 4*7^75 + 2*7^76 + 7^77 + 5*7^78 + 2*7^79 + 5*7^80 + 3*7^81 + 5*7^82 + 6*7^85 + 2*7^86 + 7^87 + 5*7^88 + 4*7^89 + 7^90 + 7^91 + 4*7^92 + 5*7^93 + 5*7^94 + 7^95 + 7^96 + 7^97 + 3*7^98 + 2*7^99 + 2*7^100 + 2*7^101 + 3*7^102 + 3*7^103 + 4*7^104 + 5*7^105 + 2*7^106 + 5*7^107 + 2*7^108 + 7^109 + 4*7^111 + 4*7^112 + 2*7^113 + 3*7^115 + 4*7^118 + 7^119 + 5*7^120 + 6*7^121 + 2*7^122 + 2*7^123 + 7^124 + 2*7^125 + 7^126 + 2*7^127 + 2*7^128 + 3*7^129 + 6*7^131 + 6*7^132 + 4*7^133 + 4*7^134 + 6*7^135 + 2*7^136 + 4*7^139 + 7^140 + 3*7^141 + 5*7^142 + 4*7^144 + 4*7^145 + 6*7^146 + 5*7^147 + 5*7^148 + 6*7^149 + O(7^150)
b = 7 + 4*7^2 + 7^3 + 2*7^4 + 5*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + 4*7^9 + 3*7^10 + 4*7^14 + 5*7^18 + 2*7^19 + 6*7^20 + 4*7^21 + 3*7^22 + 5*7^23 + 5*7^24 + 2*7^25 + 3*7^27 + 4*7^28 + 6*7^30 + 3*7^31 + 6*7^32 + 7^33 + 5*7^34 + 5*7^35 + 3*7^36 + 4*7^37 + 4*7^38 + 3*7^39 + 2*7^40 + 2*7^42 + 6*7^43 + 5*7^44 + 2*7^46 + 6*7^48 + 6*7^49 + 4*7^50 + 7^51 + 6*7^53 + 3*7^54 + 7^56 + 3*7^57 + 7^58 + 6*7^59 + 7^61 + 6*7^63 + 6*7^64 + 3*7^65 + 6*7^66 + 6*7^67 + 3*7^68 + 3*7^69 + 5*7^71 + 4*7^72 + 2*7^73 + 3*7^74 + 6*7^75 + 7^76 + 5*7^77 + 4*7^78 + 3*7^79 + 2*7^80 + 6*7^81 + 2*7^82 + 4*7^83 + 7^84 + 2*7^85 + 5*7^86 + 4*7^87 + 6*7^90 + 2*7^91 + 7^92 + 6*7^93 + 5*7^94 + 6*7^95 + 3*7^96 + 5*7^97 + 4*7^98 + 3*7^99 + 7^101 + 2*7^102 + 5*7^104 + 7^105 + 3*7^106 + 5*7^107 + 4*7^108 + 2*7^109 + 2*7^110 + 6*7^111 + 3*7^112 + 5*7^113 + 2*7^116 + 5*7^117 + 4*7^118 + 2*7^119 + 4*7^120 + 4*7^121 + 3*7^123 + 7^124 + 4*7^125 + 3*7^127 + 6*7^128 + 4*7^129 + 5*7^133 + 3*7^134 + 4*7^135 + 7^137 + 3*7^138 + 4*7^139 + 2*7^140 + 6*7^141 + 2*7^142 + 3*7^143 + 2*7^144 + 3*7^145 + 6*7^146 + 3*7^147 + 2*7^149 + O(7^150)
a = a.parent().fraction_field()(a)
b = a.parent().fraction_field()(b)
T = Matrix(ZZ,2,2,[2,2,-1,-2])
# T =Matrix(ZZ,2,2,[0,2,1,0])
# ,[j1g,j2g,j3g]

inp_vec = [(a,b,T.transpose(),qords,prec,QQ) for qords in all_possible_qords(T.transpose().change_ring(ZZ),20)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt == 'Nope':
        print 'Finished (%s) %s...'%(sum(inpt[0][3]),str(inpt[0][3]))
    else:
        try:
            i2,i4,i6,i10 = list(outt)
            assert i2 == 1 and i4 == 1 and i6 == 1
            print 'Success with %s (%s, %s, %s)'%(str(inpt[0][3]),(i2**5)/i10,(i2**3)*i4/i10,(i2**2)*i6/i10)
        except ValueError:
            print outt

# Following is some magma code:
# R<x>:=PolynomialRing(Rationals());
# C:=HyperellipticCurveFromIgusaClebsch([646, -18839/4, 4692747/8, 53905500]:Reduce:=true);
# J:=Jacobian(C);
# Cgood := HyperellipticCurve(x^6 + 6*x^5 + 11*x^4 + 14*x^3 + 5*x^2 - 12*x);
# Jgood := Jacobian(Cgood);
# EulerFactor(J,GF(31)) eq EulerFactor(Jgood,GF(31));

# IsIsomorphic(C,Cgood);
# false
# IsQuadraticTwist(C,Cgood):
# true -1



