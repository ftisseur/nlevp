function [coeffs,fun,F,sol] = qep4
%QEP4  3-by-4 QEP with known, nontrivial Jordan structure.
%  [COEFFS,FUN,F,SOL] = nlevp('qep4') generates a 3-by-4 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  The matrices are returned in a cell array: COEFFS = {A_0,A_1,A_2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  Three exact eigenvalues (0, inf, inf) and correpsonding eigenvectors
%  are returned in SOL.EVAL and  SOL.EVEC.
%  This problem has the properties pep, qep, nonregular, nonsquare,
%  real, solution.

%  Reference:
%  R. Byers, V. Mehrmann, and H. Xu. Trimmed linearizations for
%  structured matrix polynomials. Linear Algebra Appl., 429:2373-2400,
%  2008.

A_0 = [0 0 0 0; 0 0 1 0; 0 1 0 1];
A_1 = [0 1 1 0; 1 0 0 1; 1 0 0 0];
A_2 = [1 0 0 0; 0 1 0 0; 0 0 0 0];

coeffs = {A_0,A_1,A_2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
sol.eval = [0 inf inf]';
sol.evec = [2 1 0 -1; 0 0 1 0; -1 0 0 0]';

end
