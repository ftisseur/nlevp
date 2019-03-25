function [coeffs,fun,F,sol] = qep1
%QEP1  3-by-3 QEP with known eigensystem.
%  [COEFFS,FUN,F,SOL] = nlevp('qep1') generates a 3-by-3 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0
%  for which the complete eigensystem is known.
%  The matrices are returned in a cell array: COEFFS = {A_0,A_1,A_2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  The exact eigenvalues (namely 1/3, 1/2, 1, i, -i, inf) and eigenvectors
%  are known and are returned in SOL.EVAL and SOL.EVEC.
%  The eigenvalues 1/3 and 1/2 share the same eigenvector, as do i and -i.
%  This problem has the properties pep, qep, real, solution.

%  Reference:
%  F. Tisseur and K. Meerbergen, The quadratic eigenvalue problem,
%  SIAM Rev. 43 (2001), pp. 235-28 (p. 250).

A_0 = eye(3);
A_1 = [1 -6 0; 2 -7 0; 0 0 0 ];
A_2 = [0 6 0; 0 6 0 ; 0 0 1];

coeffs = {A_0,A_1,A_2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
sol.eval = [1/3 1/2 1 1i -1i inf]';
sol.evec = [1 1 0; 1 1 0; 0 1 0; 0 0 1; 0 0 1; 1 0 0]';

end
