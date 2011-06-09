# -*- coding: utf-8 -*-

"""
This TASTYL program implements the privacy-preserving set intersection protocol of:

M. J. Freedman, K. Nissim, B. Pinkas
Efficient Private Matching and Set Intersection, 
Advances in Cryptology - EUROCRYPT 2004, Vol. 3027 of LNCS, pp. 1â€“19, 2004.
http://www.pinkas.net/PAPERS/FNP04.pdf

as summarized in Sect. 3.3 of

Y. Lindell and B. Pinkas
Secure Multiparty Computation for Privacy-Preserving Data Mining 
Journal of Privacy and Confidentiality, Vol. 1, No. 1, pp. 59-98, 2009.
http://repository.cmu.edu/cgi/viewcontent.cgi?article=1004&context=jpc
"""

from tasty.crypt.math import \
    getPolyCoefficients

def protocol(c, s):
    M = 100   # size of client's set
    N = 100   # size of server's set

    c.X = ModularVec(dim=M).input(desc="X")
    s.Y = ModularVec(dim=N).input(desc="Y")

    # interpolate coeffs of poly with roots c.X
    c.a = getPolyCoefficients()(c.X)

    # encrypt and send coefficients to server
    c.ha   = HomomorphicVec(val=c.a)
    s.ha <<= c.ha

    # evaluate and rerandomize p(y_i) under enc
    s.hbarY = HomomorphicVec(dim=N)
    for i in xrange(N):   # 0, ..., N-1
        # eval poly using Horner scheme
        s.p = s.ha[M]
        for j in xrange(M-1,-1,-1):
            s.p = (s.p * s.Y[i]) + s.ha[j]
        s.hbarY[i] = s.p * Modular().rand() + \
            Homomorphic(val=s.Y[i])

    # send hbarY to client and decrypt
    c.hbarY <<= s.hbarY
    c.barY    = ModularVec(val=c.hbarY)

    # compute intersection of c.X and c.barY
    for e in c.X:
        if e in c.barY:
            c.output(e, desc="in output set")
