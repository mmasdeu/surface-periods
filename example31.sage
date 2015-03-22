load('padicperiods.sage')

# A good curve
# Conductor 31
x = QQ['x'].gen()
Crv = HyperellipticCurve((x^3 - 5*x^2 + 6*x + 1)*(x^3 - 9*x^2 +10*x - 3))
I2g, I4g, I6g, I10g = Crv.igusa_clebsch_invariants()
j1g = I2g**5 / I10g # One of absolute's Igusa invariants
j2g = I2g**3 * I4g / I10g
j3g = I2g**2 * I6g / I10g # One of absolute's Igusa invariants

p = 31
D = 1
sign = 1
prec = 40
working_prec = 80
x = QQ['x'].gen()
pol = x^2 - x - 1

path = ROOT = '/home/float/darmonpoints/'
G = BigArithGroup(p,D,1,grouptype = "SL2")
Coh = CohomologyGroup(G.Gpn)
q = ZZ(2)
K0 = (Coh.involution_at_infinity_matrix()-sign).right_kernel().change_ring(QQ)
Aq0 = Coh.hecke_matrix(q,use_magma = True).transpose().change_ring(QQ)
good_component = None
for U0,is_irred in Aq0.decomposition_of_subspace(K0):
    if U0.dimension() == 2 and is_irred and Aq0.restrict(U0).charpoly() == pol:
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
a = 11*31 + 4*31^2 + 3*31^3 + 18*31^4 + 22*31^5 + 29*31^6 + 5*31^7 + 11*31^8 + 3*31^9 + 9*31^10 + 29*31^11 + 18*31^12 + 16*31^13 + 21*31^14 + 20*31^15 + 9*31^16 + 5*31^17 + 28*31^18 + 5*31^19 + 29*31^20 + 30*31^21 + 4*31^22 + 12*31^23 + 11*31^24 + 6*31^25 + 21*31^26 + 5*31^27 + 19*31^28 + 29*31^29 + 13*31^30 + 6*31^31 + 6*31^32 + 21*31^33 + 17*31^34 + 6*31^35 + 15*31^36 + 27*31^37 + 14*31^38 + 28*31^39 + O(31^40)
b = 6*31 + 4*31^2 + 26*31^3 + 20*31^4 + 13*31^5 + 31^6 + 15*31^7 + 11*31^8 + 15*31^9 + 14*31^10 + 17*31^11 + 9*31^12 + 11*31^13 + 11*31^14 + 24*31^15 + 25*31^16 + 20*31^17 + 24*31^18 + 29*31^19 + 13*31^20 + 17*31^21 + 26*31^22 + 7*31^23 + 13*31^24 + 4*31^25 + 20*31^26 + 31^27 + 29*31^28 + 2*31^29 + 2*31^30 + 10*31^31 + 8*31^32 + 6*31^33 + 21*31^34 + 22*31^35 + 15*31^36 + 3*31^37 + 13*31^38 + 12*31^39 + O(31^40)
T = matrix(ZZ,2,2,[0,-1,-1,1])

inp_vec = [(a,b,T.transpose(),qords,30) for qords in all_possible_qords(T.transpose().change_ring(ZZ),20)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt != 'Nope':
        print 'Success with %s'%str(inpt[0][3])
        print outt
    else:
        print 'Finished %s...'%str(inpt[0][3])
