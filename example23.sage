load('padicperiods.sage')
load('darmonpoints.sage')
# A good curve
# Conductor 23
x = QQ['x'].gen()
Crv = HyperellipticCurve(x^6 - 14*x^5 + 57*x^4 -106*x^3 + 90*x^2 -16*x - 19)
I2g, I4g, I6g, I10g = Crv.igusa_clebsch_invariants()
j1g = I2g**5 / I10g # One of absolute's Igusa invariants
j2g = I2g**3 * I4g / I10g
j3g = I2g**2 * I6g / I10g # One of absolute's Igusa invariants

p = 23
D = 1
sign = 1
prec = 30
working_prec = 60
x = QQ['x'].gen()
pol = x^2 + x - 1

from sarithgroup import *
from cohomology import *

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

# Uncomment to use Teitelbaum's basis
# flist = [flist[1],-flist[0]+flist[1]]

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

# Uncomment to use Teitelbaum's basis
# g0 = ker[2]**-2*ker[3]**-2
# g1 = ker[0]

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

load('data.sage')
F = A.parent()
T = Aq0.restrict(good_component)
TF = T.change_ring(F)
a,b = p_adic_l_invariant(A,B,D,TF)
a = a.trace()/a.parent().degree()
b = b.trace()/b.parent().degree()

# Below are precomputed values
a = Qp(23,30)(7*23 + 21*23^2 + 6*23^3 + 11*23^4 + 11*23^5 + 7*23^6 + 10*23^7 + 5*23^8 + 11*23^9 + 16*23^10 + 22*23^11 + 7*23^12 + 17*23^13 + 5*23^14 + 3*23^15 + 4*23^16 + 4*23^17 + 13*23^18 + 18*23^19 + 17*23^20 + 12*23^21 + 15*23^22 + 5*23^23 + 11*23^24 + 6*23^25 + 19*23^26 + 3*23^27 + 5*23^28 + 14*23^29 + O(23^30))
b = Qp(23,30)( 18*23 + 11*23^2 + 14*23^3 + 20*23^4 + 7*23^5 + 22*23^6 + 6*23^7 + 20*23^8 + 6*23^9 + 12*23^10 + 17*23^11 + 15*23^12 + 20*23^13 + 14*23^14 + 13*23^15 + 15*23^16 + 3*23^17 + 6*23^18 + 5*23^19 + 19*23^20 + 15*23^21 + 18*23^22 + 15*23^23 + 3*23^24 + 8*23^25 + 19*23^26 + 16*23^27 + 21*23^28 + 23^29 + O(23^30))
T = matrix(ZZ,2,2,[-1,1,1,0])

inp_vec = [(a,b,T,qords,25) for qords in all_possible_qords(T.change_ring(ZZ),20)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt != 'Nope':
        print 'Success with %s'%str(inpt[0][3])
        print outt
    else:
        print 'Finished %s...'%str(inpt[0][3])
