function [coeffs,fun,F,sol] = qep3(e)
%QEP3  3-by-3 parametrized QEP with known eigensystem.
%  [COEFFS,FUN,F,SOL] = nlevp('qep3',E) generates a 3-by-3 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0
%  for which the complete eigensystem is known.
%  E is a scalar parameter, with default -1+sqrt(EPS/2);
%  The matrices are returned in a cell array: COEFFS = {A_0,A_1,A_2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  The exact eigenvalues (namely, 0, 1, 1+E, 2, 3, inf) and eigenvectors
%  are known and are returned in SOL.EVAL and SOL.EVEC.
%  This problem has the properties pep, qep, real, parameter-dependent,
%  solution.

%  Reference:
%  J.-P. Dedieu and F. Tisseur. Perturbation theory for homogeneous
%  polynomial eigenvalue problems. Linear Algebra Appl., 358:71-94, 2003
%  (p. 89).

if nargin < 1, e = -1 + sqrt(eps/2); end

A_0 = [1 -1 -1; 0 1 0; 0 0 0];
A_1 = [-3 1 0; 0 -1-e 0; 0 0 1];
A_2 = [2 0 9; 0 0 0; 0 0 -3];

coeffs = {A_2,A_1,A_0};  % Indices are reversed in the original paper.
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
sol.eval = [0 1 1+e 2 3 inf];
sol.evec = [0 1 0; 1 0 0; 1 (e-1)/(e+1) 0; 1 0 0; 0 0 1; 1 0 1]';

end
