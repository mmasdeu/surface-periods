ROOT = path = './'
# load(path + "darmonpoints.sage")
# if not path in sys.path:
#     sys.path.insert(1, path)
# del path
from itertools import product
from util import *
# from arithgroup import *
# from sarithgroup import *
# from cohomology import *

# Theta functions of (25) in Teitelbaum's
# and also the other theta functions that we need to compute the lambda's according to the formulas (23)
# We sum terms until we get the result to the accuracy prec. When deciding the set of indicies over which we sum, we have made several simplifications that probably cause that we sum more terms than strictly needed for that accuracy, but this probably won't be an issue...
def Theta(p1,p2,p3,version = None,prec = None):
    if prec is None:
        prec = 2 * p1.parent().precision_cap()
    # imax = ((1+(1+4*RR(prec/min_val)).sqrt())/2).ceiling()
    imax = (RR(1)/2 + RR(prec + RR(1)/2).sqrt()).ceiling()
    a = p1.valuation().abs()
    b = p2.valuation().abs()
    c = p3.valuation().abs()
    # print 'Entering Theta %s %s %s'%(a,b,c)
    # Define the different conditions and term forms
    if version is None:
        condition = lambda i,j: b*(i**2-ZZ(i).abs()) + a*(j**2-ZZ(j).abs()) + c*((i-j)**2 - ZZ(i-j).abs())
        resdict = {}
        resdict['1p1m'] = resdict['1p2m'] = 0
        resdict['2p2m'] = resdict['1m2p'] = 0
        resdict['3p3m'] = resdict['2m3p'] = 0
        resdict['2p3m'] = resdict['3m1p'] = resdict['3p1m'] = 0
    elif version == '1p1m':
        condition = lambda i,j: b*i**2 + a*j**2 + c*(i-j)**2
        term = lambda i,j: (-1)**j
    elif version == '1p2m':
        condition = lambda i,j:  b*(i**2-i) + a*(j**2-j) + c*(i-j)**2
        term = lambda i,j: p2**(-i)*p1**(-j)
    elif version == '2p2m':
        condition = lambda i,j: b*i**2 + a*j**2 + c*(i-j)**2
        term = lambda i,j: (-1)**i
    elif version == '1m2p':
        condition = lambda i,j: b*(i**2-i) + a*(j**2-j) + c*(i-j)**2
        term = lambda i,j: p2**(-i)*p1**(-j)*(-1)**(i+j)
    elif version == '3p3m':
        condition = lambda i,j: b*i**2 + a*j**2 + c*(i-j)**2
        term = lambda i,j: (-1)**(i+j)
    elif version == '2m3p':
        condition = lambda i,j: b*(i**2+i) + a*j**2 + c*((i-j)**2+(i-j))
        term = lambda i,j: p2**(i)*p3**((i-j))*(-1)**j
    elif version == '2p3m':
        condition = lambda i,j: b*(i**2+i) + a*j**2 + c*((i-j)**2+(i-j))
        term = lambda i,j: p2**(i)*p3**((i-j))
    elif version == '3m1p':
        condition = lambda i,j: b*i**2 + a*(j**2+j) + c*((i-j)**2+(j-i))
        term = lambda i,j: p1**(j)*p3**((j-i))*(-1)**i
    elif version == '3p1m':
        condition = lambda i,j:  b*i**2 + a*(j**2+j) + c*((i-j)**2+(j-i))
        term = lambda i,j: p1**(j)*p3**((j-i))
    else:
        raise ValueError("Wrong version? version = %s"%version)
    res = 0
    for i in range(-imax,imax + 1):
        jmax = RR(prec - .25 - i**2 +RR(i).abs())
        if jmax < 0:
            continue
        jmax = (jmax.sqrt()+.5).ceiling()
        for j in range(-jmax,jmax + 1):
            newterm = p2**(i**2) * p1**(j**2) * p3**((i-j)**2)
            if version is None:
                p1l = p1**j
                p2l = p2**i
                p3l = p3**(i-j)
                p2l_inv_p1l_inv_term = newterm * (p2l*p1l)**-1
                p2lp3l_term = newterm * p2l * p3l
                p1lp3l_inv_term = newterm * p1l * p3l**-1

                resdict['1p1m'] += newterm * (-1)**j
                resdict['1p2m'] += p2l_inv_p1l_inv_term
                resdict['2p2m'] += newterm * (-1)**i
                resdict['1m2p'] += p2l_inv_p1l_inv_term * (-1)**(i+j)
                resdict['3p3m'] += newterm * (-1)**(i+j)
                resdict['2m3p'] += p2lp3l_term * (-1)**j
                resdict['2p3m'] += p2lp3l_term
                resdict['3m1p'] += p1lp3l_inv_term * (-1)**i
                resdict['3p1m'] += p1lp3l_inv_term
            else:
                res += newterm * term(i,j)
    if version is None:
        return resdict
    else:
        return res


def lambdavec(p1, p2, p3, prec):
    th = Theta(p1,p2,p3,prec = prec)
    num = th['3p3m'] * th['2m3p']
    den = th['2p2m'] * th['2p3m']
    try:
        den.parent()._bg_ps_ring().set_default_prec(prec)
    except AttributeError:
        pass
    l1 = (num/den)**2

    num = th['1p1m'] * th['3m1p']
    den = th['3p3m'] * th['3p1m']
    try:
        den.parent()._bg_ps_ring().set_default_prec(prec)
    except AttributeError:
        pass
    l2 = (num/den)**2

    num = th['2p2m'] * th['1m2p']
    den = th['1p1m'] * th['1p2m']
    try:
        den.parent()._bg_ps_ring().set_default_prec(prec)
    except AttributeError:
        pass
    l3 = (num/den)**2
    return (l1,l2,l3)

def lambdavec_padic(p1, p2, p3):
    th = Theta(p1,p2,p3,prec = None)
    num = th['3p3m'] * th['2m3p']
    den = th['2p2m'] * th['2p3m']
    l1 = (num/den)**2
    num = th['1p1m'] * th['3m1p']
    den = th['3p3m'] * th['3p1m']
    l2 = (num/den)**2
    num = th['2p2m'] * th['1m2p']
    den = th['1p1m'] * th['1p2m']
    l3 = (num/den)**2
    return (l1,l2,l3)


def xvec(p1, p2, p3, prec):
    l1,l2,l3 = lambdavec(p1,p2,p3,prec)
    x3 = l3
    den = l2
    try:
        den.parent()._bg_ps_ring().set_default_prec(prec)
    except AttributeError:
        pass
    x2 = 1 - 1/den
    den = 1-l1
    try:
        den.parent()._bg_ps_ring().set_default_prec(prec)
    except AttributeError:
        pass
    x1 = 1/den
    return (x1,x2,x3)

def xvec_padic(p1, p2, p3):
    l1,l2,l3 = lambdavec_padic(p1,p2,p3)
    return (1/(1-l1),1 - 1/l2,l3)

def ICI_static(x1,x2,x3):
    x12, x22, x32 = x1 * x1, x2 * x2, x3 * x3
    x13, x23, x33 = x1 * x12, x2 * x22, x3 * x32
    x14, x24, x34 = x1 * x13, x2 * x23, x3 * x33
    x15, x25, x35 = x1 * x14, x2 * x24, x3 * x34
    x16, x26, x36 = x1 * x15, x2 * x25, x3 * x35
    x17, x27, x37 = x1 * x16, x2 * x26, x3 * x36
    x18, x28, x38 = x1 * x17, x2 * x27, x3 * x37

    I2 = 6*x12*x22 - 4*x12*x2*x3 - 4*x1*x22*x3 + 6*x12*x32 - 4*x1*x2*x32 + 6*x22*x32 - 4*x12*x2 - 4*x1*x22 - 4*x12*x3 + 12*x1*x2*x3 - 4*x22*x3 - 4*x1*x32 - 4*x2*x32 + 6*x12 - 4*x1*x2 + 6*x22 - 4*x1*x3 - 4*x2*x3 + 6*x32
    I4 = 4*x14*x22*x32 - 4*x13*x23*x32 + 4*x12*x24*x32 - 4*x13*x22*x33 - 4*x12*x23*x33 + 4*x12*x22*x34 - 4*x14*x22*x3 + 4*x13*x23*x3 - 4*x12*x24*x3 - 4*x14*x2*x32 + 4*x13*x22*x32 + 4*x12*x23*x32 - 4*x1*x24*x32 + 4*x13*x2*x33 + 4*x12*x22*x33 + 4*x1*x23*x33 - 4*x12*x2*x34 - 4*x1*x22*x34 + 4*x14*x22 - 4*x13*x23 + 4*x12*x24 - 4*x14*x2*x3 + 4*x13*x22*x3 + 4*x12*x23*x3 - 4*x1*x24*x3 + 4*x14*x32 + 4*x13*x2*x32 - 24*x12*x22*x32 + 4*x1*x23*x32 + 4*x24*x32 - 4*x13*x33 + 4*x12*x2*x33 + 4*x1*x22*x33 - 4*x23*x33 + 4*x12*x34 - 4*x1*x2*x34 + 4*x22*x34 - 4*x13*x22 - 4*x12*x23 + 4*x13*x2*x3 + 4*x12*x22*x3 + 4*x1*x23*x3 - 4*x13*x32 + 4*x12*x2*x32 + 4*x1*x22*x32 - 4*x23*x32 - 4*x12*x33 + 4*x1*x2*x33 - 4*x22*x33 + 4*x12*x22 - 4*x12*x2*x3 - 4*x1*x22*x3 + 4*x12*x32 - 4*x1*x2*x32 + 4*x22*x32
    I6 =8*x16*x24*x32 - 8*x15*x25*x32 + 8*x14*x26*x32 - 8*x16*x23*x33 - 4*x15*x24*x33 - 4*x14*x25*x33 - 8*x13*x26*x33 + 8*x16*x22*x34 - 4*x15*x23*x34 + 24*x14*x24*x34 - 4*x13*x25*x34 + 8*x12*x26*x34 - 8*x15*x22*x35 - 4*x14*x23*x35 - 4*x13*x24*x35 - 8*x12*x25*x35 + 8*x14*x22*x36 - 8*x13*x23*x36 + 8*x12*x24*x36 - 8*x16*x24*x3 + 8*x15*x25*x3 - 8*x14*x26*x3 - 4*x16*x23*x32 + 2*x15*x24*x32 + 2*x14*x25*x32 - 4*x13*x26*x32 - 4*x16*x22*x33 + 40*x15*x23*x33 - 28*x14*x24*x33 + 40*x13*x25*x33 - 4*x12*x26*x33 - 8*x16*x2*x34 + 2*x15*x22*x34 - 28*x14*x23*x34 - 28*x13*x24*x34 + 2*x12*x25*x34 - 8*x1*x26*x34 + 8*x15*x2*x35 + 2*x14*x22*x35 + 40*x13*x23*x35 + 2*x12*x24*x35 + 8*x1*x25*x35 - 8*x14*x2*x36 - 4*x13*x22*x36 - 4*x12*x23*x36 - 8*x1*x24*x36 + 8*x16*x24 - 8*x15*x25 + 8*x14*x26 - 4*x16*x23*x3 + 2*x15*x24*x3 + 2*x14*x25*x3 - 4*x13*x26*x3 + 24*x16*x22*x32 - 28*x15*x23*x32 + 32*x14*x24*x32 - 28*x13*x25*x32 + 24*x12*x26*x32 - 4*x16*x2*x33 - 28*x15*x22*x33 - 4*x14*x23*x33 - 4*x13*x24*x33 - 28*x12*x25*x33 - 4*x1*x26*x33 + 8*x16*x34 + 2*x15*x2*x34 + 32*x14*x22*x34 - 4*x13*x23*x34 + 32*x12*x24*x34 + 2*x1*x25*x34 + 8*x26*x34 - 8*x15*x35 + 2*x14*x2*x35 - 28*x13*x22*x35 - 28*x12*x23*x35 + 2*x1*x24*x35 - 8*x25*x35 + 8*x14*x36 - 4*x13*x2*x36 + 24*x12*x22*x36 - 4*x1*x23*x36 + 8*x24*x36 - 8*x16*x23 - 4*x15*x24 - 4*x14*x25 - 8*x13*x26 - 4*x16*x22*x3 + 40*x15*x23*x3 - 28*x14*x24*x3 + 40*x13*x25*x3 - 4*x12*x26*x3 - 4*x16*x2*x32 - 28*x15*x22*x32 - 4*x14*x23*x32 - 4*x13*x24*x32 - 28*x12*x25*x32 - 4*x1*x26*x32 - 8*x16*x33 + 40*x15*x2*x33 - 4*x14*x22*x33 + 48*x13*x23*x33 - 4*x12*x24*x33 + 40*x1*x25*x33 - 8*x26*x33 - 4*x15*x34 - 28*x14*x2*x34 - 4*x13*x22*x34 - 4*x12*x23*x34 - 28*x1*x24*x34 - 4*x25*x34 - 4*x14*x35 + 40*x13*x2*x35 - 28*x12*x22*x35 + 40*x1*x23*x35 - 4*x24*x35 - 8*x13*x36 - 4*x12*x2*x36 - 4*x1*x22*x36 - 8*x23*x36 + 8*x16*x22 - 4*x15*x23 + 24*x14*x24 - 4*x13*x25 + 8*x12*x26 - 8*x16*x2*x3 + 2*x15*x22*x3 - 28*x14*x23*x3 - 28*x13*x24*x3 + 2*x12*x25*x3 - 8*x1*x26*x3 + 8*x16*x32 + 2*x15*x2*x32 + 32*x14*x22*x32 - 4*x13*x23*x32 + 32*x12*x24*x32 + 2*x1*x25*x32 + 8*x26*x32 - 4*x15*x33 - 28*x14*x2*x33 - 4*x13*x22*x33 - 4*x12*x23*x33 - 28*x1*x24*x33 - 4*x25*x33 + 24*x14*x34 - 28*x13*x2*x34 + 32*x12*x22*x34 - 28*x1*x23*x34 + 24*x24*x34 - 4*x13*x35 + 2*x12*x2*x35 + 2*x1*x22*x35 - 4*x23*x35 + 8*x12*x36 - 8*x1*x2*x36 + 8*x22*x36 - 8*x15*x22 - 4*x14*x23 - 4*x13*x24 - 8*x12*x25 + 8*x15*x2*x3 + 2*x14*x22*x3 + 40*x13*x23*x3 + 2*x12*x24*x3 + 8*x1*x25*x3 - 8*x15*x32 + 2*x14*x2*x32 - 28*x13*x22*x32 - 28*x12*x23*x32 + 2*x1*x24*x32 - 8*x25*x32 - 4*x14*x33 + 40*x13*x2*x33 - 28*x12*x22*x33 + 40*x1*x23*x33 - 4*x24*x33 - 4*x13*x34 + 2*x12*x2*x34 + 2*x1*x22*x34 - 4*x23*x34 - 8*x12*x35 + 8*x1*x2*x35 - 8*x22*x35 + 8*x14*x22 - 8*x13*x23 + 8*x12*x24 - 8*x14*x2*x3 - 4*x13*x22*x3 - 4*x12*x23*x3 - 8*x1*x24*x3 + 8*x14*x32 - 4*x13*x2*x32 + 24*x12*x22*x32 - 4*x1*x23*x32 + 8*x24*x32 - 8*x13*x33 - 4*x12*x2*x33 - 4*x1*x22*x33 - 8*x23*x33 + 8*x12*x34 - 8*x1*x2*x34 + 8*x22*x34
    I10 = x18*x26*x34 - 2*x17*x27*x34 + x16*x28*x34 - 2*x18*x25*x35 + 2*x17*x26*x35 + 2*x16*x27*x35 - 2*x15*x28*x35 + x18*x24*x36 + 2*x17*x25*x36 - 6*x16*x26*x36 + 2*x15*x27*x36 + x14*x28*x36 - 2*x17*x24*x37 + 2*x16*x25*x37 + 2*x15*x26*x37 - 2*x14*x27*x37 + x16*x24*x38 - 2*x15*x25*x38 + x14*x26*x38 - 2*x18*x26*x33 + 4*x17*x27*x33 - 2*x16*x28*x33 + 2*x18*x25*x34 - 2*x17*x26*x34 - 2*x16*x27*x34 + 2*x15*x28*x34 + 2*x18*x24*x35 - 4*x17*x25*x35 + 4*x16*x26*x35 - 4*x15*x27*x35 + 2*x14*x28*x35 - 2*x18*x23*x36 - 2*x17*x24*x36 + 4*x16*x25*x36 + 4*x15*x26*x36 - 2*x14*x27*x36 - 2*x13*x28*x36 + 4*x17*x23*x37 - 2*x16*x24*x37 - 4*x15*x25*x37 - 2*x14*x26*x37 + 4*x13*x27*x37 - 2*x16*x23*x38 + 2*x15*x24*x38 + 2*x14*x25*x38 - 2*x13*x26*x38 + x18*x26*x32 - 2*x17*x27*x32 + x16*x28*x32 + 2*x18*x25*x33 - 2*x17*x26*x33 - 2*x16*x27*x33 + 2*x15*x28*x33 - 6*x18*x24*x34 + 4*x17*x25*x34 + 4*x16*x26*x34 + 4*x15*x27*x34 - 6*x14*x28*x34 + 2*x18*x23*x35 + 4*x17*x24*x35 - 6*x16*x25*x35 - 6*x15*x26*x35 + 4*x14*x27*x35 + 2*x13*x28*x35 + x18*x22*x36 - 2*x17*x23*x36 + 4*x16*x24*x36 - 6*x15*x25*x36 + 4*x14*x26*x36 - 2*x13*x27*x36 + x12*x28*x36 - 2*x17*x22*x37 - 2*x16*x23*x37 + 4*x15*x24*x37 + 4*x14*x25*x37 - 2*x13*x26*x37 - 2*x12*x27*x37 + x16*x22*x38 + 2*x15*x23*x38 - 6*x14*x24*x38 + 2*x13*x25*x38 + x12*x26*x38 - 2*x18*x25*x32 + 2*x17*x26*x32 + 2*x16*x27*x32 - 2*x15*x28*x32 + 2*x18*x24*x33 - 4*x17*x25*x33 + 4*x16*x26*x33 - 4*x15*x27*x33 + 2*x14*x28*x33 + 2*x18*x23*x34 + 4*x17*x24*x34 - 6*x16*x25*x34 - 6*x15*x26*x34 + 4*x14*x27*x34 + 2*x13*x28*x34 - 2*x18*x22*x35 - 4*x17*x23*x35 - 6*x16*x24*x35 + 24*x15*x25*x35 - 6*x14*x26*x35 - 4*x13*x27*x35 - 2*x12*x28*x35 + 2*x17*x22*x36 + 4*x16*x23*x36 - 6*x15*x24*x36 - 6*x14*x25*x36 + 4*x13*x26*x36 + 2*x12*x27*x36 + 2*x16*x22*x37 - 4*x15*x23*x37 + 4*x14*x24*x37 - 4*x13*x25*x37 + 2*x12*x26*x37 - 2*x15*x22*x38 + 2*x14*x23*x38 + 2*x13*x24*x38 - 2*x12*x25*x38 + x18*x24*x32 + 2*x17*x25*x32 - 6*x16*x26*x32 + 2*x15*x27*x32 + x14*x28*x32 - 2*x18*x23*x33 - 2*x17*x24*x33 + 4*x16*x25*x33 + 4*x15*x26*x33 - 2*x14*x27*x33 - 2*x13*x28*x33 + x18*x22*x34 - 2*x17*x23*x34 + 4*x16*x24*x34 - 6*x15*x25*x34 + 4*x14*x26*x34 - 2*x13*x27*x34 + x12*x28*x34 + 2*x17*x22*x35 + 4*x16*x23*x35 - 6*x15*x24*x35 - 6*x14*x25*x35 + 4*x13*x26*x35 + 2*x12*x27*x35 - 6*x16*x22*x36 + 4*x15*x23*x36 + 4*x14*x24*x36 + 4*x13*x25*x36 - 6*x12*x26*x36 + 2*x15*x22*x37 - 2*x14*x23*x37 - 2*x13*x24*x37 + 2*x12*x25*x37 + x14*x22*x38 - 2*x13*x23*x38 + x12*x24*x38 - 2*x17*x24*x32 + 2*x16*x25*x32 + 2*x15*x26*x32 - 2*x14*x27*x32 + 4*x17*x23*x33 - 2*x16*x24*x33 - 4*x15*x25*x33 - 2*x14*x26*x33 + 4*x13*x27*x33 - 2*x17*x22*x34 - 2*x16*x23*x34 + 4*x15*x24*x34 + 4*x14*x25*x34 - 2*x13*x26*x34 - 2*x12*x27*x34 + 2*x16*x22*x35 - 4*x15*x23*x35 + 4*x14*x24*x35 - 4*x13*x25*x35 + 2*x12*x26*x35 + 2*x15*x22*x36 - 2*x14*x23*x36 - 2*x13*x24*x36 + 2*x12*x25*x36 - 2*x14*x22*x37 + 4*x13*x23*x37 - 2*x12*x24*x37 + x16*x24*x32 - 2*x15*x25*x32 + x14*x26*x32 - 2*x16*x23*x33 + 2*x15*x24*x33 + 2*x14*x25*x33 - 2*x13*x26*x33 + x16*x22*x34 + 2*x15*x23*x34 - 6*x14*x24*x34 + 2*x13*x25*x34 + x12*x26*x34 - 2*x15*x22*x35 + 2*x14*x23*x35 + 2*x13*x24*x35 - 2*x12*x25*x35 + x14*x22*x36 - 2*x13*x23*x36 + x12*x24*x36
    return I2, I4, I6, I10

def IgusaClebschFromHalfPeriods(a, b, c, prec = None):
    if a.valuation() == 0:
        a, c = c, a
    elif b.valuation() == 0:
        b, c = c, b
    if a.valuation() < 0:
        a,b,c = 1/a, 1/b, 1/c
    if a.valuation() <= 0 or b.valuation() <= 0 or c.valuation() < 0:
        raise RuntimeError
    if prec is None:
        return ICI_static(*xvec_padic(a,b,c))
    else:
        return ICI_static(*xvec(a,b,c,prec))

# computes the p-adic L-invariant
# A = <gamma_1,gamma_1>
# B = <gamma_1,gamma_2>
# D = <gamma_2,gamma_2>
# Tmatrix is the matrix of the T_ell operator with respect to the basis (gamma_1,gamma_2)
# the output is (a,b), where L_p = a + bT
def p_adic_l_invariant(A,B,D,Tmatrix):
    K = A.parent()
    x,y,z,t = Tmatrix.change_ring(K).list()
    alpha,beta,delta = A.ordp(),B.ordp(),D.ordp()
    M = Matrix(K,3,2,[alpha,x*alpha+z*beta,beta,y*alpha+t*beta,delta,y*beta+t*delta])
    n  = Matrix(K,3,1,[A.log(0),B.log(0),D.log(0)])
    return M.solve_right(n).list()

def qlogs_from_Lp_and_ords(a,b,Tmatrix,q1ord, q2ord, q3ord):
    K = a.parent()
    ordA = q2ord + q3ord
    ordB = -q3ord
    ordD = q1ord + q3ord
    bord = matrix(K,2,2,[ordA,ordB,ordB,ordD]) * Tmatrix
    M = Matrix(K,3,2,[ordA, bord[0,0], ordB, bord[0,1],ordD,bord[1,1]])
    logA, logB, logD = (M * matrix(K,2,1,[a,b])).list()
    return logB + logD, logA + logB, -logB

def all_possible_qords(Tmatrix,N,initial = None):
    # Find linear equation from matrix
    x, y, z, t = [ZZ(o) for o in Tmatrix.list()]
    r = x+y-z-t
    ans = set([])
    if z != 0:
        # Know that q1^z = q2^y q3^r
        for ordq2, ordq3 in product(range(N+1),repeat = 2):
            ordq1 = ZZ(y * ordq2 + r * ordq3)
            if ordq1 % z != 0:
                continue
            ordq1 /= z
            ordtup = [ordq1,ordq2,ordq3]
            if all([o >= 0 for o in ordtup]) and len([o for o in ordtup if o == 0]) <= 1:
                ans = ans.union(ans,set([(ordq1, ordq2, ordq3)]))
    if y != 0:
        # Know that q2^y = q1^z q3^-r
        for ordq1, ordq3 in product(range(N+1),repeat = 2):
            ordq2 = ZZ(z * ordq1 - r * ordq3)
            if ordq2 % y != 0:
                continue
            ordq2 /= y
            ordtup = [ordq1,ordq2,ordq3]
            if all([o >= 0 for o in ordtup]) and len([o for o in ordtup if o == 0]) <= 1:
                ans = ans.union(ans,set([(ordq1, ordq2, ordq3)]))
    if r != 0:
        # Know that q3^r = q1^z q2^-y
        for ordq1, ordq2 in product(range(N+1),repeat = 2):
            ordq3 = ZZ(z * ordq1 - y * ordq2)
            if ordq3 % r != 0:
                continue
            ordq3 /= r
            ordtup = [ordq1,ordq2,ordq3]
            if all([o >= 0 for o in ordtup]) and len([o for o in ordtup if o == 0]) <= 1:
                ans = ans.union(ans,set([(ordq1, ordq2, ordq3)]))
    tmp = sorted(list(ans),key = lambda x: x[0]+x[1]+x[2])
    t0 = 0
    if initial is not None:
        try:
            t0 = tmp.index(initial)
        except ValueError:
            t0 = 0
    return tmp[t0:]

# use hensel lemma to lift an approximate root x0 of the polynomial f to a root to the desired precision
def ComputeRootFromApprox(f,x0, prec):
    xn = x0
    while True:
        xnn = xn - f.subs(xn)/f.derivative().subs(xn)
        if (xnn-xn).valuation() > prec:
            break
        xn = xnn
    return xn

def recognize_j1(j1,base = QQ,phi = None,threshold = 0.9):
    deg = base.degree()
    Kp = j1.parent()
    p = Kp.prime()
    if phi is None:
        phi = lambda t:t
    threshold = threshold * RR(p).log(10)
    j1 = p**-j1.valuation() * j1
    fx = algdep(j1,deg)
    if height_polynomial(fx) < threshold * j1.precision_relative():
        return fx
    raise ValueError('Unrecognized')

def recognize_invariants(j1,j2,j3,pval,base = QQ,phi = None):
    deg = base.degree()
    x = QQ['x'].gen()
    Kp = j1.parent()
    j2, j3 = Kp(j2), Kp(j3)
    p = Kp.prime()
    if phi is None:
        phi = lambda t:t
    threshold = .9 * RR(p).log(10)
    for froot,_ in algdep(p**pval * j1,deg).roots(base):
        if froot.global_height() > threshold * j1.precision_relative():
            continue
        # At this point we may have recognized j1
        j1 = p**-pval * froot
        I10 = p**pval * froot.denominator()
        verbose('I10 = %s'%I10)
        I2_5 = froot.denominator() * froot
        try:
            for q,e in I2_5.factor():
                if e % 5 != 0:
                    v = 5 - (e % 5)
                    qv = q**v
                    I2_5 *= qv
                    I10 *= qv
        except ArithmeticError: pass
        verbose('I2_5 = %s, I10 = %s'%(I2_5,I10))
        for I2,_ in (x**5 - I2_5).roots(base):
            verbose('I2 = %s'%I2)
            I4p = phi(I10/I2**3) * j2
            for I4,_ in algdep(I4p,deg).roots(base):
                verbose('I4 = %s'%I4)
                if I4.global_height() > threshold * I4p.precision_relative():
                    continue
                I6p = phi(I10/I2**2) * j3
                for I6,_ in algdep(I6p,deg).roots(base):
                    verbose('I6 = %s'%I6)
                    if I6.global_height() > threshold * I6p.precision_relative():
                        continue
                    return (I2, I4, I6, I10)
    raise ValueError('Unrecognized')

@parallel
def find_igusa_invariants_from_L_inv(a,b,T,qords,prec,base = QQ,cheatjs = None,phi = None):
    F = a.parent()
    TF = T.change_ring(F)
    p = F.prime()
    x = QQ['x'].gen()
    K.<y> = F.extension(x^2 + p)
    deg = base.degree()
    oq1, oq2, oq3 = qords
    q10, q20, q30 = [o.exp() for o in qlogs_from_Lp_and_ords(a,b,TF,oq1,oq2,oq3)]
    for s1, s2, s3 in product(F.teichmuller_system(),repeat = 3):
        q1 = K(s1 * q10 * p**oq1)
        q2 = K(s2 * q20 * p**oq2)
        q3 = K(s3 * q30 * p**oq3)
        for p1,p2,p3 in product(our_sqrt(q1,K,return_all = True),our_sqrt(q2,K,return_all = True),our_sqrt(q3,K,return_all = True)):
            try:
                # p1,p2,p3 = our_sqrt(q1,K),our_sqrt(q2,K),our_sqrt(q3,K)
                I2c, I4c, I6c, I10c = IgusaClebschFromHalfPeriods(p1,p2,p3,prec = prec)
                # Get absolute invariants j1, j2, j3
                j1 = I2c**5/I10c
                j2 = I2c**3*I4c/I10c
                j3 = I2c**2*I6c/I10c
                j1n = j1.trace() / j1.parent().degree()
                j2n = j2.trace() / j2.parent().degree()
                j3n = j3.trace() / j3.parent().degree()
                assert (j1 - j1n).valuation() > 5,'j1 = %s, j1n = %s'%(j1,j1n)
                assert (j2 - j2n).valuation() > 5,'j2err'
                assert (j3 - j3n).valuation() > 5,'j3err'
                j1, j2, j3 = j1n, j2n, j3n

                if cheatjs is not None:
                    if min([(u-v).valuation() for u,v in zip([j1,j2,j3],cheatjs)]) > 3:
                        return (oq1,oq2,oq3,1)
                else:
                    # return recognize_invariants(j1,j2,j3,oq1+oq2+oq3,base = base,phi = phi)
                    return (1,1,1,recognize_j1(j1,base = base,phi = phi,threshold = 0.90))
            except ValueError:
                pass
            except Exception as e:
                return str(e.message)
    return 'Nope'

def euler_factor_twodim(p,T):
    x = QQ['x'].gen()
    t = T.trace()
    n = T.determinant()
    return x**4 - t*x**3 + (2*p+n)*x**2 - p*t*x + p*p
