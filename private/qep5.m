function [coeffs,fun,F] = qep5
%QEP5  3-by-3 nonregular QEP with known Smith form.
%  [COEFFS,FUN,F] = nlevp('qep5') generates a 3-by-3 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  Q has Smith form diag(1,lambda-1,0).
%  The matrices are returned in a cell array: COEFFS = {A_0,A_1,A_2}.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, nonregular, real.

%  Reference:
%  Paul M. Van Dooren and Patrick Dewilde, The eigenstructure of an
%  arbitrary polynomial matrix: computational aspects, Linear Algebra
%  Appl., 50:545-579, 1983.

A_0 = [1 2 -2; 0 -1 -2; 0 0 0];
A_1 = [1 3 0; 1 4 2; 0 -1 -2];
A_2 = [1 4 2; 0 0 0; 1 4 2];

coeffs = {A_0,A_1,A_2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
