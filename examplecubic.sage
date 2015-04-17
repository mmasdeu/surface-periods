load('padicperiods.sage')

sign = 1
prec = 30
working_prec = 2000

path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

#### Search for candidates
# candidates = load('all_candidates_with_covol.sobj')
# for ii,cand in enumerate(candidates):
#     print 'Doing ii = %s'%ii
#     if ii < 69:
#         continue
#     F,P,D,Np,Sinf,covol = cand
#     if F.degree() == 2:
#         continue
#     Sinf_places = [v for v,o in zip(F.real_places(prec = Infinity),Sinf) if o == -1]
#     try:
#         abtuple = quaternion_algebra_invariants_from_ramification(F,D,Sinf_places)
#     except ValueError:
#         continue
#     G = BigArithGroup(p,abtuple,Np,base = F,grouptype = 'PGL2')
#     Coh = CohomologyGroup(G.Gpn)
#     if Coh.dimension() >= 2:
#         fwrite(str((ii,F.polynomial(),p.gens_reduced()[0],D.gens_reduced()[0],Np.gens_reduced()[0],Sinf)),'g2candidates.txt')

# PARAMETERS

x = QQ['x'].gen()
F.<r> = NumberField(x^3 - x**2 + 2*x - 3)
P = F.ideal(r)
D = F.ideal(r^2 - r + 1)
Np = F.ideal(1)
Sinf = [-1]
Sinf_places = [v for v,o in zip(F.real_places(prec = Infinity),Sinf) if o == -1]
abtuple = quaternion_algebra_invariants_from_ramification(F,D,Sinf_places)

G = BigArithGroup(P,abtuple,Np,base = F,grouptype = 'PGL2')
Coh = CohomologyGroup(G.Gpn)

flist, hecke_data = Coh.get_twodim_cocycle(sign,return_all = False)
ell, T = hecke_data[0]
g0, g1 = G.get_pseudo_orthonormal_homology(flist,smoothen = ell)

from homology import *
xi10,xi20 = lattice_homology_cycle(G,g0,working_prec)
xi11,xi21 = lattice_homology_cycle(G,g1,working_prec)

Phif = get_overconvergent_class_quaternionic(p,flist[0],G,prec,sign,progress_bar = True)
Phig = get_overconvergent_class_quaternionic(p,flist[1],G,prec,sign,progress_bar = True)

from integrals import *
A = integrate_H1_many(G,[(xi10,False,1),(xi20,True,-1)],Phif,1,method = 'moments',prec = working_prec, progress_bar = True)
B = integrate_H1_many(G,[(xi11,False,1),(xi21,True,-1)],Phif,1,method = 'moments',prec = working_prec, progress_bar = True)
D = integrate_H1_many(G,[(xi11,False,1),(xi21,True,-1)],Phig,1,method = 'moments',prec = working_prec, progress_bar = True)

A = (A + O(p**(prec+A.valuation()))).trace()/A.parent().degree()
B = (B + O(p**(prec+B.valuation()))).trace()/B.parent().degree()
D = (D + O(p**(prec+D.valuation()))).trace()/D.parent().degree()

F = A.parent()
TF = T.change_ring(F)
a1,b1 = p_adic_l_invariant(A,B,D,TF)

# g0 = ker[3]**2
# g1 = ker[0]**2*ker[3]
# print [flist[0].evaluate(g0),flist[0].evaluate(g1)]
# print [flist[1].evaluate(g0),flist[1].evaluate(g1)]


# Below are precomputed values
a = 3 + 3^3 + 3^4 + 2*3^7 + 3^9 + 3^10 + 3^11 + 2*3^12 + 3^14 + 3^15 + 2*3^17 + 3^18 + 2*3^20 + 3^22 + 3^23 + 3^24 + 3^27 + 2*3^28 + 2*3^29 + 3^30 + 3^31 + 2*3^32 + 2*3^36 + 2*3^38 + 3^39 + 3^40 + 3^41 + 2*3^45 + 3^46 + 3^47 + 2*3^48 + 2*3^50 + 3^52 + 2*3^55 + 3^56 + 3^57 + 2*3^58 + 2*3^59 + 3^60 + 3^64 + 2*3^69 + 3^70 + 2*3^71 + 2*3^73 + 2*3^74 + 3^76 + 3^77 + 3^79 + 3^80 + 3^81 + 2*3^82 + 3^83 + 3^85 + 2*3^86 + 2*3^87 + 2*3^88 + 2*3^89 + 2*3^90 + 2*3^91 + 2*3^92 + 3^95 + 2*3^98 + 2*3^100 + 2*3^102 + 2*3^103 + 3^106 + 2*3^107 + 2*3^108 + 2*3^109 + 2*3^110 + 3^111 + 2*3^112 + 2*3^114 + 2*3^116 + 3^117 + 2*3^119 + 3^120 + 3^121 + 3^123 + 2*3^124 + 3^126 + 3^127 + 2*3^129 + 2*3^132 + 3^133 + 3^134 + 3^135 + 2*3^137 + 3^138 + 2*3^139 + 3^141 + 3^142 + 2*3^144 + 2*3^145 + 2*3^146 + 3^148 + 3^150 + 2*3^151 + 2*3^153 + 3^154 + 2*3^155 + 3^156 + 3^158 + 3^159 + 2*3^160 + 2*3^161 + 2*3^162 + 2*3^163 + 3^165 + 2*3^166 + 2*3^167 + 3^170 + 3^172 + 3^174 + 2*3^175 + 3^176 + 2*3^177 + 3^178 + 3^180 + 3^181 + 2*3^182 + 2*3^183 + 2*3^184 + 2*3^185 + 3^186 + 3^187 + 3^188 + 2*3^189 + 2*3^190 + 2*3^191 + 3^192 + 3^193 + O(3^195)
b = 2*3 + 3^5 + 3^6 + 3^7 + 3^8 + 3^9 + 3^10 + 2*3^11 + 3^12 + 3^13 + 3^14 + 3^15 + 3^16 + 3^18 + 2*3^19 + 3^20 + 3^21 + 3^22 + 3^23 + 3^25 + 2*3^26 + 2*3^29 + 2*3^30 + 3^34 + 3^35 + 3^36 + 3^37 + 3^38 + 2*3^39 + 2*3^41 + 2*3^42 + 3^44 + 3^46 + 3^47 + 2*3^48 + 3^49 + 3^51 + 2*3^54 + 2*3^55 + 3^56 + 3^58 + 2*3^60 + 2*3^65 + 2*3^68 + 3^69 + 2*3^71 + 2*3^73 + 2*3^76 + 3^79 + 2*3^80 + 3^82 + 3^83 + 3^84 + 2*3^85 + 2*3^87 + 2*3^88 + 3^89 + 2*3^94 + 3^96 + 3^97 + 2*3^98 + 3^100 + 2*3^104 + 2*3^106 + 2*3^107 + 2*3^108 + 3^109 + 3^110 + 2*3^111 + 2*3^112 + 3^114 + 2*3^117 + 3^118 + 2*3^119 + 3^120 + 3^121 + 3^122 + 2*3^124 + 3^125 + 3^128 + 3^129 + 2*3^130 + 3^131 + 3^135 + 3^136 + 2*3^137 + 2*3^139 + 2*3^140 + 3^141 + 2*3^142 + 2*3^143 + 3^145 + 2*3^146 + 2*3^150 + 3^151 + 2*3^152 + 3^153 + 3^155 + 2*3^156 + 3^157 + 2*3^160 + 2*3^161 + 2*3^163 + 2*3^167 + 2*3^170 + 3^173 + 2*3^175 + 2*3^176 + 3^177 + 3^178 + 3^182 + 3^184 + 3^185 + 2*3^186 + 3^187 + 2*3^188 + 3^190 + 2*3^192 + 3^194 + O(3^196)
T =Matrix(ZZ,2,2,[-3,-5,-2,-4])

inp_vec = [(a,b,T.transpose(),qords,195,P.ring()) for qords in all_possible_qords(T.transpose().change_ring(ZZ),10)]
for inpt in inp_vec:
    find_igusa_invariants_from_L_inv(*inpt)

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



