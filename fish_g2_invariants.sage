
load_attach_path('/home/float/GitProjects/surface-periods')
load('padicperiods.sage')
sign = 1
prec = 20
timeout = 12 * 3600 # 12 hour timeout

#path = ROOT = '/home/float/darmonpoints/'
from sarithgroup import *
from cohomology import *

x = QQ['x'].gen()
r = QQ['r'].gen()
# PARAMETERS
# This file will try to compute many p-adic L-invariants, hoping that some work.
load('candidates_atr.sage')
print 'Candidates LOADED'

input_vec = []
for pol, Pgen, Dgen, Npgen, Pnm, Dnm, dim in data:
    input_vec.append((0,pol,Pgen,Dgen,Npgen,[-1 for o in range(pol.degree() - 1)]))

f_guess_equation = fork(guess_equation,timeout = timeout)
for inpt, outp in parallel(lambda a,b,c,d,e,f:f_guess_equation(a,b,c,d,e,f,sign,prec))(input_vec):
    print 'Finished inpt = %s. Result: %s'%(str(inpt),str(outp))
