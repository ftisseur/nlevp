function [coeffs,fun,F] = time_delay
%TIME_DELAY    3-by-3 NEP from time-delay system.
%   [COEFFS,FUN,F] = nlevp('time_delay') generates a 3-by-3
%   nonlinear matrix function -lambda*EYE(3) + A_0 + A_1*exp(-lambda)
%   that is the characteristic equation of a time-delay system with a
%   single delay and constant coefficients.
%   The problem has a double non-semisimple eigenvalue lambda = 3*PI*i.
%   The matrices are returned in a cell array: COEFFS = {EYE(3), A_0, A_1}.
%   FUN is a function handle to evaluate the functions -lambda, 1, and
%   exp(-lambda) and their derivatives.
%   F is the function handle -lambda*EYE(3) + A_0 + A_1*exp(-lambda).
%   This problem has the properties nep, real.

%   A typo in coefficient b(2) has been corrected as in Jarlebring and
%   Michiels (2011).

%   Reference:
%   E. Jarlebring. Convergence factors of {Newton} methods for nonlinear
%      eigenvalue problems. Linear Algebra Appl., 2011. In press, corrected
%      proof. DOI:10.1016/j.laa.2010.08.045.
%   E. Jarlebring and W. Michiels. Invariance properties in the root
%      sensitivity of time-delay systems with double imaginary roots,
%      Automatica: 46 (2010), pp. 1112-1115.
%   E. Jarlebring and W. Michiels.
%      Analyzing the convergence factor of residual inverse iteration.
%      BIT, 2011,  In press, corrected proof. DOI:10.1007/s10543-011-0336-2.

a(1) = 2*(65*pi+32)/(5*(8+5*pi));
a(2) = 9*pi^2*(13+5*pi)/(8+5*pi);
a(3) = 324*pi^2*(5*pi+4)/(5*(8+5*pi));

b(1) = (260*pi + 128 + 225*pi^2)/ (10*(8+5*pi));
b(2) = 45*pi^2/(8+5*pi);  % Typo corrected from references above.
b(3) = 81*pi^2 *(40*pi + 32 + 25*pi^2)/ (10*(8+5*pi));

A_0 = [0 1 0; 0 0 1; -a(3) -a(2) -a(1)];
A_1 = [0 0 0; 0 0 0; -b(3) -b(2) -b(1)];

coeffs = {eye(3), A_0, A_1};
fun = @(lam) time_delay_fun(lam);
F = @(lam) -lam*coeffs{1}+coeffs{2}+coeffs{3}*exp(-lam);

end


function varargout = time_delay_fun(lam)

lam = lam(:);
n = length(lam);

varargout{1} = [-lam, ones(n,1), exp(-lam)];
if nargout >= 2
    f = -exp(-lam);
    varargout{2} = [-ones(n,1),zeros(n,1),f];
    for i=2:nargout-1
        f = -f;
        varargout{i+1} = [zeros(n,2),f];
    end
end


end