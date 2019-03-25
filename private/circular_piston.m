function [coeffs,fun,F] = circular_piston
%CIRCULAR_PISTON  Sparse QEP from model of circular piston.
%  [COEFFS,FUN,F] = nlevp('circular_piston') generates the coefficient
%  matrices of a quadratic matrix polynomial lambda^2*M + lambda*E + K
%  from an axi-symmetric infinite element model for a circular piston.
%  This problem has dimension 2025. It is nonsymmetric with a singular
%  leading coefficient matrix M.
%  The matrices are sparse and returned in a cell array: COEFFS = {K, E, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*E + lambda^2*M.
%  This problem has the properties pep, qep, real, sparse.

%  Reference:
%  J.-P. Coyette, K. Meerbergen, and M. Robb. Time integration for spherical 
%  acoustic finite-infinite element models. Internat. J. Numer. Methods Eng., 
%  64(13):1752–1768, 2005.

load('piston.mat')

coeffs = {K,E,M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
