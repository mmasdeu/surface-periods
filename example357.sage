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
q = ZZ(5)
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
a = 4*7 + 7^2 + 5*7^3 + 4*7^4 + 7^5 + 4*7^7 + 2*7^8 + 7^11 + 4*7^13 + 2*7^14 + 5*7^15 + 6*7^16 + 2*7^17 + 4*7^18 + 4*7^19 + 5*7^20 + 7^21 + 4*7^22 + 6*7^23 + 2*7^25 + 6*7^26 + 6*7^27 + 2*7^28 + 7^31 + 4*7^32 + 3*7^33 + 3*7^34 + 7^35 + 6*7^37 + 7^38 + 3*7^39 + 5*7^40 + 7^41 + 2*7^42 + 7^43 + 3*7^44 + 7^45 + 7^46 + 3*7^47 + 7^49 + 6*7^50 + 5*7^51 + 5*7^52 + 6*7^53 + 5*7^55 + 3*7^56 + 5*7^57 + 4*7^58 + 3*7^59 + 3*7^60 + 4*7^61 + 3*7^62 + 7^63 + 3*7^64 + 4*7^65 + 2*7^67 + 5*7^68 + 4*7^69 + O(7^70)
b = 7 + 4*7^2 + 7^3 + 2*7^4 + 5*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + 4*7^9 + 3*7^10 + 4*7^14 + 5*7^18 + 2*7^19 + 6*7^20 + 4*7^21 + 3*7^22 + 5*7^23 + 5*7^24 + 2*7^25 + 3*7^27 + 4*7^28 + 6*7^30 + 3*7^31 + 6*7^32 + 7^33 + 5*7^34 + 5*7^35 + 3*7^36 + 4*7^37 + 4*7^38 + 3*7^39 + 2*7^40 + 2*7^42 + 6*7^43 + 5*7^44 + 2*7^46 + 6*7^48 + 6*7^49 + 4*7^50 + 7^51 + 6*7^53 + 3*7^54 + 7^56 + 3*7^57 + 7^58 + 6*7^59 + 7^61 + 6*7^63 + 6*7^64 + 3*7^65 + 6*7^66 + 6*7^67 + 3*7^68 + 3*7^69 + O(7^70)
T =Matrix(QQ,2,2,[0,2,1,0])

inp_vec = [(a,b,T.transpose(),qords,prec-5,QQ,[j1g,j2g,j3g]) for qords in all_possible_qords(T.transpose().change_ring(ZZ),20)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt != 'Nope':
        i2,i4,i6,i10 = list(outt)
        print 'Success with %s (%s, %s, %s)'%(str(inpt[0][3]),i2**5/i10,i2**3*i4/i10,i2**2*i6/i10)
        print outt
    else:
        print 'Finished %s...'%str(inpt[0][3])

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



