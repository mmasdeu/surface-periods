load('padicperiods.sage')

# A good curve
# Conductor 165 = 3 * 5 * 11
x = QQ['x'].gen()
Crv = HyperellipticCurve(x^6 + 6*x^5 + 11*x^4 + 14*x^3 + 5*x^2 - 12*x)
I2g, I4g, I6g, I10g = Crv.igusa_clebsch_invariants()
j1g = I2g**5 / I10g # One of absolute's Igusa invariants
j2g = I2g**3 * I4g / I10g
j3g = I2g**2 * I6g / I10g # One of absolute's Igusa invariants

p = 11
D = 1
Np = 15
sign = 1
prec = 10
working_prec = 120
x = QQ['x'].gen()
pol = x^2 + 2*x - 1

path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *
G = BigArithGroup(p,D,Np)
Coh = CohomologyGroup(G.Gpn)


def get_pseudo_orthonormal_homology(self, cocycles, smoothen = 0):
    ker = self.get_homology_kernel(smoothen = smoothen)
    # Now we need to find elements orthogonal to given cocycles
    dim = len(cocycles)
    A = Matrix(ZZ,dim,0)
    for vec0 in ker:
        A = A.augment(vector([ZZ(f.evaluate(vec0)[0]) for f in cocycles]))
    B = A.solve_right(Matrix(ZZ,dim,dim,1))
    B = B.denominator() * B
    Gab = self.Gpn.abelianization()
    return [ Gab.ab_to_G(sum(ZZ(i) * Gab.G_to_ab(self.Gpn(o)) for o,i in zip(ker,col.list()))).quaternion_rep for col in B.columns() ]

flist, hecke_data = Coh.get_twodim_cocycle(sign,pol = pol,return_all = False)
ell, T = hecke_data[0]
g0, g1 = get_pseudo_orthonormal_homology(G,flist,smoothen = ZZ(ell))

from homology import *
xi10,xi20 = lattice_homology_cycle(G,g0,working_prec)
xi11,xi21 = lattice_homology_cycle(G,g1,working_prec)
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
a = 3*11 + 8*11^2 + 3*11^4 + 9*11^5 + 7*11^6 + 4*11^7 + 9*11^8 + 9*11^9 + 2*11^11 + 4*11^12 + 11^13 + 3*11^14 + 3*11^15 + 3*11^16 + 10*11^17 + 3*11^18 + 3*11^19 + 3*11^20 + 3*11^21 + 11^22 + 11^23 + 2*11^24 + 4*11^25 + 9*11^26 + 7*11^27 + 9*11^28 + 5*11^29 + 11^30 + 7*11^31 + 7*11^32 + 4*11^34 + 5*11^35 + 3*11^36 + 4*11^37 + 3*11^38 + 4*11^39 + 8*11^40 + 4*11^41 + 11^42 + 3*11^43 + 2*11^44 + 10*11^45 + 4*11^47 + 3*11^48 + 3*11^49 + 5*11^50 + 2*11^51 + 10*11^52 + 9*11^53 + 7*11^54 + 8*11^55 + 9*11^56 + 2*11^58 + 6*11^59 + O(11^60)
b = 7*11^2 + 3*11^3 + 5*11^4 + 2*11^5 + 4*11^6 + 8*11^8 + 11^9 + 8*11^10 + 7*11^11 + 2*11^12 + 4*11^13 + 9*11^14 + 11^15 + 11^16 + 5*11^18 + 10*11^19 + 6*11^21 + 6*11^22 + 10*11^24 + 9*11^25 + 3*11^26 + 11^27 + 6*11^28 + 9*11^29 + 7*11^30 + 10*11^31 + 5*11^32 + 11^33 + 6*11^34 + 7*11^35 + 8*11^37 + 11^38 + 10*11^39 + 9*11^40 + 8*11^41 + 2*11^42 + 3*11^43 + 10*11^44 + 10*11^46 + 4*11^47 + 2*11^49 + 7*11^50 + 10*11^51 + 4*11^54 + 2*11^55 + 8*11^56 + 7*11^57 + 4*11^58 + O(11^60)
T =Matrix(ZZ,2,2,[-3,-2,1,1])


inp_vec = [(a,b,T.transpose(),qords,50) for qords in all_possible_qords(T.transpose().change_ring(ZZ),10)]
for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
    if outt != 'Nope':
        i2,i4,i6,i10 = list(outt)
        print 'Success with %s (%s, %s, %s)'%(str(inpt[0][3]),i2**5/i10,i2**3*i4/i10,i2**2*i6/i10)
        print outt
    else:
        print 'Finished %s...'%str(inpt[0][3])


# Here is the result:
# There are more successes, this is the first (not the lowest height one)
# Success with (0, 2, 2)
# (16490, 1161961/4, 12864116757/8, 39530700) # Not the right curve
# Success with (1, 0, 2)
# (646, -18839/4, 4692747/8, 53905500)

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



