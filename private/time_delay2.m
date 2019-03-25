function [coeffs,fun,F] = time_delay2(tau)
%TIME_DELAY2  2-by-2 NEP from a time-delay system.
%  [COEFFS,FUN,F] = nlevp('time_delay2', tau) generates a 2-by-2
%  nonlinear matrix function lambda*EYE(2) + B_0 + A_1*exp(-tau*lambda)
%  that is the characteristic equation of the time-delay system x'(t) =
%  B_0*x(t) + A_1*x(t-tau).
%  The matrices are returned in a cell array: COEFFS = {EYE(2), B_0, A_1}.
%  FUN is a function handle to evaluate the functions lambda, 1, and
%  exp(tau*lambda) and their derivatives.
%  F is the function handle lambda*EYE(2) + B_0 + A_1*exp(tau*lambda).
%  This problem has the properties nep, parameter-dependent, real.

%  Reference:
%  W. Michiels and S.-I. Niculescu. Stability and stabilization of time-delay 
%  systems. vol. 12 of Advances in Design and Control, Society for Industrial 
%  and Applied Mathematics (SIAM), Philadelphia, PA, 2007.
%
%  D. Kressner, A block Newton method for nonlinear eigenvalue problems,
%  Numer. Math., 114:355-372, 2009.

if nargin < 1 || isempty(tau)
tau = 1;
end

B_0 = [5 -1; -2 6];
A_1 = [2 -1; -4 1];

coeffs = {eye(2), B_0, A_1};
fun = @(lam) time_delay2_fun(tau, lam);
F = @(lam) lam*coeffs{1} + coeffs{2}+coeffs{3}*exp(-tau*lam);

end


function varargout = time_delay2_fun(tau, lam)

lam = lam(:);
n = length(lam);

varargout{1} = [lam, ones(n,1), exp(-tau*lam)];
if nargout >= 2
    f = -tau*exp(-tau*lam);
    varargout{2} = [ones(n,1), zeros(n,1), f];
    for i=2:nargout-1
        f = -tau*f;
        varargout{i+1} = [zeros(n,2),f];
    end
end
end
