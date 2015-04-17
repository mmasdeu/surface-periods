load('padicperiods.sage')

sign = 1
prec = 30
working_prec = 3000

path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

# Search for candidates
candidates = load('all_candidates_with_covol.sobj')
for ii,cand in enumerate(candidates):
    print 'Doing ii = %s'%ii
    if ii < 29:
        continue
    F,P,D,Np,Sinf,covol = cand
    if F.degree() == 2:
        continue
    Sinf_places = [v for v,o in zip(F.real_places(prec = Infinity),Sinf) if o == -1]
    try:
        abtuple = quaternion_algebra_invariants_from_ramification(F,D,Sinf_places)
    except ValueError:
        continue
    G = BigArithGroup(P,abtuple,Np,base = F,grouptype = 'PGL2')
    Coh = CohomologyGroup(G.Gpn)
    if Coh.dimension() >= 2:
        print 'Success (%s: %s)'%(ii,F.polynomial())
        fwrite(str((ii,F.polynomial(),P.gens_reduced()[0],D.gens_reduced()[0],Np.gens_reduced()[0],Sinf)),'g2candidates.txt')
