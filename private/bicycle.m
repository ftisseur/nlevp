function [coeffs,fun,F] = bicycle(v)
%BICYCLE  2-by-2 QEP from the Whipple bicycle model.
%  [COEFFS,FUN,F] = nlevp('bicycle',v) constructs a 2-by-2 quadratic matrix
%  polynomial lambda^2*M + lambda*v*C + K0 + v^2*K2 arising in the
%  study of the dynamic behaviour of a bicycle.
%  v (default 5) is the forward speed in m/s.
%  The matrices are returned in a cell array: COEFFS = {K, C, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*C + lambda^2*M.
%  This problem has the properties pep, qep, real, parameter_dependent.

%  Reference:
%  J. P. Meijaard, J. M. Papadopoulos, A. Ruina and A. L. Schwab, Linearized
%  dynamics equations for the balance and steer of a bycicle: a benchmark
%  and review, Proc. Roy. Soc. London Ser. A, 463(2084):1955-1982, 2007.

if nargin < 1 || isempty(v), v = 5; end
g = 9.81; % N/kg

M = [80.81722 2.31941332208709; 2.31941332208709 0.29784188199686];
K0 = [-80.95 -2.59951685249872; -2.59951685249872 -0.80329488458618];
K2 = [0  76.59734589573222; 0 2.65431523794604];
C = [0 33.86641391492494; -0.85035641456978 1.68540397397560];

coeffs = {g*K0+v^2*K2, v*C, M};

fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
