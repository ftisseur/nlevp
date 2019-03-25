function [coeffs,fun,F] = mobile_manipulator
%MOBILE_MANIPULATOR  QEP from model of 2-dimensional 3-link mobile manipulator.
%  [COEFFS,FUN,F] = nlevp('mobile_manipulator') generates a 5x5 quadratic
%  matrix polynomial lambda^2*M + lambda*D + K arising from the modelling
%  as a time-invariant descriptor control system of a two-dimensional
%  three-link mobile manipulator.
%  The matrices are returned in a cell array: coeffs = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*M + lambda*D + K.
%  This problem has the properties pep, qep, real.

%  References:
%  R. Byers, C. He, and V. Mehrmann, Where is the Nearest
%     Non-Regular Pencil?, Linear Algebra Appl., 285 (1998), pp. 81-105.
%  N. J. Higham, and F. Tisseur, More on pseudospectra for polynomial
%     eigenvalue problems and applications in control theory,
%     Linear Algebra Appl., 351-352 (2002), pp. 435-453.

K0 = [67.4894 69.2393 -69.2393
      69.8124 1.68624 -1.68617
     -69.8123 -1.68617 -68.2707];

M0 = [18.7532 -7.94493 7.94494
     -7.94493  31.8182 -26.8182
      7.94494 -26.8182  26.8182];

D0 = [-1.52143 -1.55168  1.55168
       3.22064  3.28467 -3.28467
      -3.22064 -3.28467  3.28467];

F0 = [1 0 0
      0 0 1];

coeffs{3} = [M0 zeros(3,2); zeros(2,5)];
coeffs{2} = [D0 zeros(3,2); zeros(2,5)];
coeffs{1} = [K0 -F0'; F0 zeros(2,2)];

fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
