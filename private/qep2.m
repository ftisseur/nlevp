function [coeffs,fun,F,sol] = qep2
%QEP2  3-by-3 QEP with known, nontrivial Jordan structure.
%  [COEFFS,FUN,F,SOL] = nlevp('qep2') generates a 3-by-3 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0
%  for which the complete eigensystem is known.
%  The matrices are returned in a cell array: COEFFS = {A_0,A_1,A_2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  The exact eigenvalues (namely -1, 1, 1, 1, inf, inf) and correspomnding
%  Jordan chains are known and are returned in SOL.EVAL and SOL.EVEC.
%  This problem has the properties pep, qep, real, solution.

%  Reference:
%  F. Tisseur and K. Meerbergen, The quadratic eigenvalue problem,
%  SIAM Rev. 43 (2001), pp. 235-286 (p. 256).

A_0 = diag([1 -1 1]);
A_1 = [-2 0 1; 0 0 0; 0 0 0];
A_2 = [1 0 0; 0 1 0; 0 0 1];

coeffs = {A_0,A_1,A_2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
sol.eval = [-1 1 1 1 inf inf]';
sol.evec = [0 1 0; 0 1 0; 1 0 0; 0 1 0; 0 0 1; -1 0 1]';

end
