function [coeffs,fun,F] = closed_loop(alpha)
%CLOSED_LOOP  2-by-2 QEP associated with closed-loop control system.
%  [COEFFS,FUN,F] = nlevp('closed_loop',ALPHA) constructs a 2-by-2 quadratic
%  matrix polynomial lambda^2*A2 + lambda*A1 + A0 associated with a
%  closed-loop control system with feedback gains 1 and 1 + ALPHA.
%  For 0 < ALPHA < 0.875 all eigenvalues lie inside the unit disc.
%  The default is ALPHA = 1. The matrices are returned in a cell array:
%  COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle A_0 + lambda*A_1 + lambda^2*A_2.
%  This problem has the properties pep, qep, real, parameter_dependent.

%  Reference:
%  Francoise Tisseur and Nicholas J. Higham. Structured pseudospectra for
%  polynomial eigenvalue problems, with applications. SIAM J. Matrix
%  Anal. Appl., 23(1):187-208, 2001.

if nargin < 1 || isempty(alpha), alpha = 1; end

coeffs{3} = eye(2);
coeffs{2} = [0 1+alpha; 1 0];
coeffs{1} = diag([0.5 0.25]);

fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
