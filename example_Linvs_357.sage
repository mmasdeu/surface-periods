load('padicperiods.sage')
path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

# set_verbose(1)
# Set a particular seed
# magma.eval("SetSeed(1251311619)")
# mseed = ZZ.random_element(10000)
mseed = 8033
magma.eval("SetSeed(%s)"%mseed)
print 'mseed = %s'%mseed

prec = 10
working_prec = 400
sign = 1
x = QQ['x'].gen()
pol = x^2 - 2

# First calculation
p = 7
D = 1
Np = 3*17

G = BigArithGroup(p,D,Np,grouptype = 'PSL2')
Coh = CohomologyGroup(G.Gpn)
flist, hecke_data = Coh.get_twodim_cocycle(sign,pol = pol,return_all = False)
ell, T = hecke_data[0]
g0, g1 = G.get_pseudo_orthonormal_homology(flist,smoothen = ell)

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
TF = T.change_ring(F)
a1,b1 = p_adic_l_invariant(A,B,D,TF)

## Second calculation
p = 7
D = 3*17
Np = 1

G = BigArithGroup(p,D,Np,grouptype = 'PSL2')
Coh = CohomologyGroup(G.Gpn)
flist, hecke_data = Coh.get_twodim_cocycle(sign,pol = pol,return_all = False)
ell, T = hecke_data[0]


def get_pseudo_orthonormal_homology(self, cocycles, smoothen,ker,A):

    dim = len(cocycles)
    A = Matrix(ZZ,dim,0)
    for vec0 in ker:
        A = A.augment(vector([ZZ(f.evaluate(vec0)[0]) for f in cocycles]))
    Gab = self.Gpn.abelianization()
    kernrms = [ Gab.G_to_ab(self.Gpn(o)).vector() for o in ker]
    def custom_norm(B):
        Gab = self.Gpn.abelianization()
        return max([(sum(ZZ(i) * o for o,i in zip(kernrms,col.list()))).norm(1) for col in B.columns()])

    min_norm = 10000
    minB = None
    for i in range(len(A.columns())):
        B0 = A.submatrix(col = i,ncols = 1)
        if B0.is_zero():
            continue
        for j in range(i+1,len(A.columns())):
            if A.column(j).is_zero():
                continue
            B1 = B0.augment(A.column(j))
            if B1.determinant().is_zero():
                continue
            B = B1.adjoint()
            B = B.parent()(1/B.gcd() * B)
            Bnrm = custom_norm(B)
            print Bnrm
            if Bnrm < min_norm:
                min_norm = Bnrm
                minB = B
                minij = (i,j)
    assert minB is not None
    B  = minB
    ans = []
    i,j = minij
    ker = [ker[i],ker[j]]
    Gab = self.Gpn.abelianization()
    return [ Gab.ab_to_G(sum(ZZ(i) * Gab.G_to_ab(self.Gpn(o)) for o,i in zip(ker,col.list()))).quaternion_rep for col in B.columns() ]

g0, g1 = get_pseudo_orthonormal_homology(G,flist,ZZ(ell),ker,A)

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
TF = T.change_ring(F)
a2,b2 = p_adic_l_invariant(A,B,D,TF)

print 'RESULTS'
print '-------'
print 'a1 = %s'%a1
print 'a2 = %s'%a2
print 'Error = p^%s'%(a1-a2).valuation()
print 'b1 = %s'%b1
print 'b2 = %s'%b2
print 'Error = p^%s'%(b1-b2).valuation()

# Final output
# RESULTS
# -------
# a1 = 4*7 + 7^2 + 5*7^3 + 4*7^4 + 7^5 + 4*7^7 + 2*7^8 + O(7^9)
# a2 = 4*7 + 7^2 + 5*7^3 + 4*7^4 + 7^5 + 4*7^7 + 2*7^8 + O(7^9)
# Error = p^9
# b1 = 7 + 4*7^2 + 7^3 + 2*7^4 + 5*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + O(7^9)
# b2 = 7 + 4*7^2 + 7^3 + 2*7^4 + 5*7^5 + 6*7^6 + 6*7^7 + 6*7^8 + O(7^9)
# Error = p^9
