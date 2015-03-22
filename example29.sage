load('padicperiods.sage')

# A good curve
# Conductor 29
x = QQ['x'].gen()
Crv = HyperellipticCurve(x^6 - 4*x^5 - 12*x^4 + 2*x^3 + 8*x^2 + 8*x - 7)
I2g, I4g, I6g, I10g = Crv.igusa_clebsch_invariants()
j1g = I2g**5 / I10g # One of absolute's Igusa invariants
j2g = I2g**3 * I4g / I10g
j3g = I2g**2 * I6g / I10g # One of absolute's Igusa invariants

p = 29
D = 1
sign = 1
prec = 40
working_prec = 80
x = QQ['x'].gen()
pol = x^2 + 2*x - 1

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
a = 8*29 + 23*29^2 + 2*29^3 + 21*29^4 + 17*29^5 + 13*29^6 + 26*29^7 + 8*29^8 + 17*29^9 + 16*29^10 + 13*29^11 + 24*29^12 + 18*29^13 + 26*29^14 + 11*29^15 + 2*29^16 + 29^17 + 2*29^18 + 8*29^19 + 11*29^20 + 28*29^22 + 2*29^25 + 24*29^26 + 14*29^27 + 13*29^28 + 3*29^29 + 5*29^30 + 23*29^31 + 6*29^32 + 21*29^33 + 20*29^34 + 20*29^35 + 5*29^36 + 22*29^37 + 17*29^38 + 24*29^39 + O(29^40)
b = 12*29 + 3*29^3 + 12*29^4 + 6*29^5 + 26*29^6 + 6*29^7 + 15*29^8 + 7*29^9 + 28*29^10 + 7*29^11 + 8*29^12 + 15*29^13 + 2*29^14 + 4*29^15 + 15*29^16 + 24*29^17 + 18*29^18 + 28*29^19 + 23*29^20 + 27*29^21 + 22*29^22 + 2*29^23 + 12*29^24 + 9*29^25 + 27*29^26 + 18*29^27 + 21*29^28 + 13*29^29 + 20*29^30 + 5*29^31 + 13*29^32 + 14*29^33 + 15*29^34 + 22*29^35 + 26*29^36 + 23*29^37 + 11*29^38 + 5*29^39 + O(29^40)
T = matrix(ZZ,2,2,[-2,1,1,0])
# T = matrix(ZZ,2,2,[1,-2,1,-3]) # Teitelbaum
inp_vec = [(a,b,T,qords,30) for qords in all_possible_qords(T.change_ring(ZZ),20)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt != 'Nope':
        print 'Success with %s'%str(inpt[0][3])
        print outt
    else:
        print 'Finished %s...'%str(inpt[0][3])
