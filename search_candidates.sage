load('padicperiods.sage')

sign = 1
prec = 30
working_prec = 3000

path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

# Search for candidates
candidates = load('all_candidates_with_covol.sobj')

@fork(timeout = 20 * 60)
def dimension_test(P,abtuple,Np,F):
    G = BigArithGroup(P,abtuple,Np,base = F,grouptype = 'PGL2')
    Coh = CohomologyGroup(G.Gpn)
    return Coh.dimension()

@parallel
def try_candidate(ii,cand):
    # print 'Doing ii = %s'%ii
    F,P,D,Np,Sinf,covol = cand
    if F.degree() == 2:
        return
    Sinf_places = [v for v,o in zip(F.real_places(prec = Infinity),Sinf) if o == -1]
    try:
        abtuple = quaternion_algebra_invariants_from_ramification(F,D,Sinf_places)
    except ValueError:
        return
    try:
        dim = ZZ(dimension_test(P,abtuple,Np,F))
    except TypeError:
        return
    if dim >= 2:
        # print 'Success (%s: %s)'%(ii,F.polynomial())
        fwrite(str((ii,F.polynomial(),P.gens_reduced()[0],D.gens_reduced()[0],Np.gens_reduced()[0],Sinf)),'g2candidates.txt')
