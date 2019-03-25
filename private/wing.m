function [coeffs,fun,F] = wing
%WING     3-by-3 QEP from analysis of oscillations of a wing in an airstream.
%  [COEFFS,FUN] = nlevp('wing') generates a 3-by-3 quadratic matrix
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0 arising from
%  the analysis of the oscillations of a wing in an airstream.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2,
%  and their derivatives.
%  F is the function handleQ(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  This problem has the properties pep, qep, real.

%  References:
%  R. A. Frazer, W. J. Duncan, and A. R. Collar, Elementary Matrices
%     and Some Applications to Dynamics and Differential Equations,
%     Cambridge University Press, tenth ed., 1938.  1963 printing.
%  P. Lancaster, Lambda-Matrices and Vibrating Systems, Pergamon Press,
%     Oxford, 1966, xiii+196 pp. Reprinted by Dover, New York, 2002.
%  F. Tisseur and N. J. Higham, Structured Pseudospectra for Polynomial
%     Eigenvalue Problems, with Applications, SIAM J. Matrix Anal. Appl.,
%     23(1):187-208, 2001.

coeffs{3} = [17.6 1.28  2.89;
             1.28 0.824 0.413;
             2.89 0.413 0.725];

coeffs{2} = [7.66 2.45  2.1;
             0.23 1.04  0.223;
             0.6  0.756 0.658];

coeffs{1} = [121  18.9 15.9;
             0    2.7   0.145;
             11.9 3.64 15.5];

fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end