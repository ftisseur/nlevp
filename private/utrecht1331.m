function [coeffs,fun,F] = utrecht1331
%UTRECHT1331  1331-by-1331 QEP from propagation of sound waves.
%  [COEFFS,FUN,F] = nlevp('utrecht1331') generates a 1331-by-1331 quadratic matrix
%  polynomial Q(lambda) = K + lambda*D + lambda^2*M, where K and M
%  are real symmetric positive definite, while D is complex,
%  non-Hermitian and singular.
%  The matrices are returned in a cell array: COEFFS = {K,D,M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = K + lambda*D + lambda^2*M.
%  This problem has the properties pep, qep, sparse, banded (132 bands).

%  Reference:
%  G. Sleijpen, H. van der Vorst and M. van Gijzen, Quadratic
%  Eigenproblems are No Problems, SIAM News, 8:9–10, September 1996.


load utrecht1331

coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
