# -*- coding: utf-8 -*-

"""
This TASTYL program implements the efficient privacy-preserving
face recognition protocol based on the Eigenface algorithm.

For details on the protocol, see the original paper:

Ahmad-Reza Sadeghi, Thomas Schneider, and Immo Wehrenberg.
Efficient privacy-preserving face recognition.
In 12th International Conference on Information Security and Cryptology (ICISC'09),
Volume 5984 of LNCS. Springer, December 2-4, 2009.
Full version available at http://eprint.iacr.org/2009/507.
"""

def protocol(c, s):
  K = 12      # dimension of eigenspace
  N = 10304   # number of pixels
  M = 42      # size of database

  # Declarations
  s.homegabar = HomomorphicVec(dim=K)
  s.hgamma    = HomomorphicVec(dim=N)
  s.hD        = HomomorphicVec(dim=M)
  c.bot = Unsigned(val=M,bitlen=bitlength(M+1))
  c.gbot = Garbled(val=c.bot)

  # Client inputs
  c.gamma=UnsignedVec(bitlen=8, dim=N).input()

  # Server inputs
  s.omega = UnsignedVec(bitlen=32, dim=(M, K)).input()
  s.psi = UnsignedVec(bitlen=8, dim=N).input()
  s.u = SignedVec(bitlen=8, dim=(K, N)).input()
  s.tau = Unsigned(bitlen=50).input()

  # Projection
  s.hgamma <<= HomomorphicVec(val=c.gamma)
  for i in xrange(K):
    s.homegabar[i] = Homomorphic(val=-(s.u[i].\
      dot(s.psi)))+ (s.hgamma.dot(s.u[i]))

  # Distance
  s.hs3 = s.homegabar.dot(s.homegabar)
  for i in xrange(M):
    s.hD[i] = s.hs3 + s.omega[i].dot(s.omega[i])
    s.hD[i] += s.homegabar.dot(s.omega[i]*(-2))

  # Minimum
  c.gD <<= GarbledVec(val=s.hD, force_bitlen=50, force_signed=False)
  c.gDmin_val,c.gDmin_ix=c.gD.min_value_index()
  c.gtau <<= Garbled(val=s.tau)
  c.gcmp = c.gDmin_val < c.gtau
  c.gout = c.gcmp.mux(c.gbot, c.gDmin_ix)
  c.out = Unsigned(val=c.gout)
  if c.out == c.bot:
    c.output("no match found")
  else:
    c.out.output(desc="matched index in DB")
