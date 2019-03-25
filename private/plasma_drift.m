function [coeffs,fun,F] = plasma_drift(n)
%PLASMA_DRIFT Cubic PEP arising in Tokamak reactor design.
%  [COEFFS,FUN,F] = nlevp('plasma_drift',N) generates the coefficient matrices
%  of a cubic matrix polynomial lambda^3*M3 + lambda^2*M2 + lambda*M1
%  +  M0 of dimension N modelling drift instabilities in the plasma
%  edge inside a Tokamak reactor.  The physically interesting eigenvalue
%  is the one with largest imaginary part.
%  Only two values of N are allowed: N = 128 (the default) and N = 512.
%  The matrices are returned in a cell array: COEFFS = {M0, M1, M2, M3}.
%  FUN is a function handle to evaluate the monomials 1, lambda, lambda^2,
%  lambda^3 and their derivatives.
%  F is the function handle lambda^3*M3 + lambda^2*M2 + lambda*M1 + M0.
%  This problem has the property pep and, if N = 512, sparse.

%  Reference:
%  M. Tokar, F. Kelly and X. Loozen, Role of thermal instabilities and
%  anomalous transport in threshold of detachment and mulitfacetted
%  asymmetric radiation from the edge (MARFE), Physics of Plasmas 12,
%  052510 (2005),

if nargin < 1, n = 128; end

if n == 128
   load plasma_drift_128
  elseif n == 512
   load plasma_drift_512
else
   error('This value of N not supported.')
end

coeffs = {M0,M1,M2,M3};
fun = @(lam) nlevp_monomials(lam,3);
F = @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*(coeffs{3}+lam*coeffs{4}));

end
