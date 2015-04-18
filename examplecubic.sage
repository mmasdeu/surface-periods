load('padicperiods.sage')

sign = 1
prec = 300
working_prec = 3000

#path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

x = QQ['x'].gen()
r = QQ['r'].gen()
# PARAMETERS
candidate_list =[('129_0', x^3 - x^2 + 2*x - 3, r, r^2 - r + 1, 1, [-1]),
                 ('129_1', x^3 - x^2 + 2*x - 3, r^2 - r + 1, r, 1, [-1]),
                 ('132_0', x^3 - x^2 + 4*x - 1, r^2 - r + 3, -r^2 + r - 2, 1, [-1]),
                 ('132_1', x^3 - x^2 + 4*x - 1, -r^2 + r - 2, r^2 - r + 3, 1, [-1]),
                 (130, x^3 - x^2 + 2*x - 3, r, 2, 1, [-1]),
                 (183, x^3 - x^2 - 3*x - 3, -r^2 + 2*r + 2, r + 1, 1, [-1]),
                 (184, x^3 - x^2 - 3*x - 3, -r^2 + 2*r + 2, r, 1, [-1]),
                 (207, x^3 - x^2 + 5*x - 2, r, -r + 1, 1, [-1]),
                 (194, x^3 - x^2 - 2*x - 3, r, -r + 1, 1, [-1]),
                 (209, x^3 - x^2 + 5*x - 2, 2*r^2 - r + 9, -r + 1, 1, [-1]),
                 (210, x^3 - 4*x - 5, -r - 1, -r^2 + r + 3, 1, [-1]),
                 (236, x^3 - x^2 - 5*x + 8, -r + 2, -r^2 - r + 3, 1, [-1]),
                 (269, x^3 - x^2 - 3*x - 4, r^2 - 2*r - 2, -r - 1, 1, [-1])]
@parallel
def try_your_luck(code,pol,Pgen,Dgen,Npgen,Sinf):
    outfile = 'out_try_luck_cubic_%s_%s.txt'%(i,code)
    fwrite('Starting computation for candidate %s'%str((code,pol,Pgen,Dgen,Npgen,Sinf)),outfile)
    F.<r> = NumberField(pol)
    r = F.gen()
    P = F.ideal(Pgen)
    D = F.ideal(Dgen)
    Np = F.ideal(Npgen)
    Sinf_places = [v for v,o in zip(F.real_places(prec = Infinity),Sinf) if o == -1]
    abtuple = quaternion_algebra_invariants_from_ramification(F,D,Sinf_places)

    G = BigArithGroup(P,abtuple,Np,base = F,grouptype = 'PGL2')
    Coh = CohomologyGroup(G.Gpn)
    fwrite('Computed Cohomology group',outfile)
    flist, hecke_data = Coh.get_twodim_cocycle(sign,return_all = False)
    fwrite('Obtained cocycle',outfile)
    ell, T = hecke_data[0]
    g0, g1 = G.get_pseudo_orthonormal_homology(flist,smoothen = ell)

    fwrite('Obtained homology generators',outfile)
    from homology import *
    xi10,xi20 = lattice_homology_cycle(G,g0,working_prec)
    xi11,xi21 = lattice_homology_cycle(G,g1,working_prec)
    fwrite('Defined homology cycles',outfile)
    Phif = get_overconvergent_class_quaternionic(P,flist[0],G,prec,sign,progress_bar = True)
    Phig = get_overconvergent_class_quaternionic(P,flist[1],G,prec,sign,progress_bar = True)
    fwrite('Overconvergent lift completed',outfile)
    from integrals import *
    num = integrate_H1(G,xi10,Phif,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
    den = integrate_H1(G,xi20,Phif,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
    A = num/den
    fwrite('Finished computation of A period',outfile)
    num = integrate_H1(G,xi11,Phif,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
    den = integrate_H1(G,xi21,Phif,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
    B = num/den
    fwrite('Finished computation of B period',outfile)

    num = integrate_H1(G,xi11,Phig,1,method = 'moments',prec = working_prec, twist = False,progress_bar = True)
    den = integrate_H1(G,xi21,Phig,1,method = 'moments',prec = working_prec, twist = True,progress_bar = True)
    D = num/den
    fwrite('Finished computation of D period',outfile)

    A = A.add_bigoh(prec + A.valuation())
    B = B.add_bigoh(prec + B.valuation())
    D = D.add_bigoh(prec + D.valuation())

    A = A.trace()/A.parent().degree()
    B = B.trace()/B.parent().degree()
    D = D.trace()/D.parent().degree()

    fwrite('A = %s'%A,outfile)
    fwrite('B = %s'%A,outfile)
    fwrite('D = %s'%A,outfile)
    fwrite('T = %s'%str(T.list()),outfile)
    F = A.parent()
    TF = T.change_ring(F)
    a,b = p_adic_l_invariant(A,B,D,TF)

    fwrite('a = %s'%a,outfile)
    fwrite('b = %s'%b,outfile)

    fwrite('Trying to recognize invariants...',outfile)
    phi = G._F_to_local
    inp_vec = [(a,b,T.transpose(),qords,prec,P.ring(),None,phi) for qords in all_possible_qords(T.transpose().change_ring(ZZ),20)]
    for inpt in inp_vec:
        ans = find_igusa_invariants_from_L_inv(*inpt)
        if ans != 'Nope':
            fwrite(str(ans),outfile)
    fwrite('Done',outfile)


for inpt, outp in try_your_luck(candidate_list):
    print 'Finished inpt = %s'%str(inpt)

# # Below are precomputed values
# a = 4*7 + 2*7^2 + 5*7^3 + 3*7^4 + 3*7^5 + 4*7^6 + 2*7^7 + 3*7^8 + 6*7^9 + 2*7^10 + 4*7^11 + 3*7^13 + 5*7^14 + 4*7^15 + 3*7^16 + 3*7^17 + 2*7^18 + 2*7^19 + 6*7^20 + 3*7^21 + 7^22 + 2*7^23 + 5*7^24 + 4*7^25 + 4*7^26 + 7^28 + 6*7^29 + 6*7^30 + 7^31 + 5*7^32 + 6*7^33 + 3*7^35 + 5*7^36 + 6*7^37 + 6*7^38 + 3*7^39 + 7^40 + 3*7^41 + 4*7^42 + 2*7^43 + 3*7^44 + 2*7^45 + 2*7^46 + 2*7^47 + 7^48 + 7^49 + 5*7^50 + 2*7^51 + 6*7^52 + 5*7^53 + 7^54 + 3*7^55 + 4*7^57 + 7^58 + 6*7^59 + 2*7^60 + 5*7^62 + 3*7^63 + 4*7^64 + 3*7^65 + 4*7^66 + 5*7^68 + 2*7^69 + 5*7^70 + 5*7^71 + 7^72 + 4*7^73 + 7^74 + 6*7^75 + 4*7^76 + 2*7^77 + 5*7^78 + 7^79 + 4*7^80 + 6*7^81 + 2*7^82 + 4*7^83 + 7^84 + 3*7^85 + 5*7^86 + 4*7^87 + 4*7^89 + 6*7^90 + 4*7^91 + 6*7^92 + 7^93 + 4*7^94 + 3*7^95 + 4*7^96 + 7^99 + 5*7^101 + 6*7^102 + 4*7^104 + 2*7^105 + 7^106 + 5*7^107 + 3*7^108 + 6*7^109 + 7^110 + 6*7^111 + 6*7^112 + 4*7^113 + 6*7^114 + 5*7^115 + 2*7^116 + 2*7^118 + 5*7^119 + 6*7^120 + 5*7^123 + 6*7^124 + 3*7^125 + 5*7^127 + 4*7^128 + 3*7^131 + 5*7^132 + 5*7^133 + 5*7^134 + 4*7^135 + 5*7^137 + 4*7^138 + 6*7^139 + 5*7^140 + 2*7^141 + 3*7^143 + 2*7^144 + 5*7^145 + 2*7^146 + 5*7^147 + 7^148 + 2*7^149 + 7^150 + 4*7^151 + 4*7^152 + 3*7^154 + 3*7^155 + 7^156 + 2*7^157 + 5*7^158 + 3*7^159 + 2*7^160 + 2*7^161 + 4*7^162 + 4*7^163 + 3*7^164 + 6*7^167 + 3*7^168 + 2*7^169 + 4*7^171 + 7^172 + 7^173 + 5*7^174 + 4*7^175 + 2*7^176 + 6*7^177 + 7^178 + 4*7^179 + 6*7^180 + 6*7^181 + 5*7^182 + 4*7^184 + 6*7^185 + 2*7^186 + 3*7^188 + 4*7^189 + 4*7^190 + 5*7^191 + 7^192 + 3*7^193 + 3*7^195 + 6*7^196 + 3*7^197 + O(7^200)
# b = 7^2 + 6*7^4 + 2*7^5 + 2*7^6 + 6*7^7 + 5*7^8 + 2*7^9 + 6*7^10 + 3*7^11 + 2*7^12 + 3*7^13 + 2*7^14 + 2*7^15 + 5*7^16 + 2*7^18 + 5*7^19 + 5*7^20 + 5*7^21 + 7^22 + 7^24 + 2*7^25 + 3*7^26 + 6*7^27 + 7^28 + 7^29 + 5*7^30 + 3*7^32 + 6*7^33 + 6*7^34 + 7^35 + 4*7^36 + 5*7^37 + 6*7^38 + 5*7^39 + 2*7^40 + 6*7^41 + 2*7^42 + 6*7^45 + 2*7^46 + 2*7^47 + 2*7^48 + 5*7^49 + 5*7^50 + 3*7^51 + 2*7^52 + 6*7^53 + 2*7^54 + 5*7^55 + 3*7^56 + 5*7^57 + 5*7^58 + 5*7^59 + 5*7^60 + 5*7^61 + 6*7^62 + 6*7^63 + 4*7^64 + 2*7^65 + 7^66 + 4*7^67 + 4*7^68 + 2*7^69 + 6*7^70 + 7^71 + 6*7^72 + 6*7^73 + 7^75 + 2*7^76 + 6*7^77 + 3*7^78 + 2*7^79 + 7^82 + 5*7^83 + 5*7^84 + 6*7^85 + 2*7^86 + 2*7^87 + 3*7^88 + 4*7^89 + 3*7^90 + 3*7^92 + 7^93 + 2*7^94 + 4*7^95 + 5*7^96 + 6*7^97 + 7^98 + 4*7^99 + 6*7^100 + 6*7^101 + 3*7^103 + 5*7^104 + 5*7^105 + 5*7^106 + 2*7^107 + 3*7^109 + 5*7^111 + 4*7^112 + 4*7^113 + 4*7^115 + 7^116 + 6*7^117 + 2*7^118 + 2*7^120 + 5*7^122 + 2*7^123 + 7^124 + 2*7^125 + 7^126 + 7^127 + 2*7^128 + 2*7^129 + 5*7^130 + 6*7^131 + 3*7^132 + 6*7^134 + 3*7^135 + 5*7^136 + 5*7^137 + 2*7^138 + 4*7^139 + 4*7^141 + 3*7^142 + 6*7^143 + 6*7^144 + 3*7^147 + 6*7^148 + 4*7^150 + 4*7^151 + 6*7^152 + 7^153 + 6*7^154 + 6*7^156 + 7^157 + 6*7^158 + 6*7^159 + 2*7^160 + 5*7^161 + 5*7^162 + 3*7^163 + 5*7^164 + 2*7^165 + 2*7^166 + 3*7^167 + 2*7^168 + 2*7^169 + 7^170 + 2*7^171 + 5*7^172 + 5*7^173 + 7^174 + 4*7^175 + 2*7^176 + 5*7^178 + 2*7^179 + 2*7^180 + 2*7^181 + 3*7^182 + 4*7^183 + 6*7^184 + 2*7^185 + 4*7^186 + 6*7^187 + 3*7^188 + 4*7^189 + 5*7^190 + 2*7^191 + 7^192 + 7^193 + 2*7^194 + 5*7^195 + 5*7^196 + 7^197 + 2*7^199 + O(7^200)
# T =Matrix(ZZ,2,2,[-1, 2, 5, 0])

# for inpt, outt in find_igusa_invariants_from_L_inv(inp_vec):
#     if outt != 'Nope':
#         try:
#             i2,i4,i6,i10 = list(outt)
#             print 'Success with %s (%s, %s, %s)'%(str(inpt[0][3]),i2**5/i10,i2**3*i4/i10,i2**2*i6/i10)
#         except ValueError:
#             print outt
#     else:
#         print 'Finished %s...'%str(inpt[0][3])


