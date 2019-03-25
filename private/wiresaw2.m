function [coeffs,fun,F] = wiresaw2(n,v,eta)
%WIRESAW2  QEP from vibration analysis of wiresaw with viscous damping.
%  [COEFFS,FUN,F] = nlevp('wiresaw2',N,V,ETA) constructs an N-by-N quadratic
%  polynomial Q(lambda) = lambda^2*M + lambda*D + K
%  arising in the vibration analysis of wiresaw with viscous damping effects
%  on the wire. This is the problem from NLEVP('WIRESAW1',N,V) with
%  added viscous damping.
%  V and ETA are nonnegative parameters corresponding to the
%  speed of the wire and damping, repectively.
%  The default values are N = 10, V = 0.01 and ETA = 0.8.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2,
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*M + lambda*D + K.
%  This problem has the properties pep, qep, real, parameter-dependent,
%  scalable.

%  Reference:
%  S. Wei and I. Kao, Vibration analysis of wire and frequency response in
%  the modern wiresaw manufacturing process, Journal of Sound and Vibration
%  213(5) (2000), pp. 1383-1395.

if nargin < 3 || isempty(eta), eta = 0.8; end;
if nargin < 2 || isempty(v), v = 0.01; end;
if nargin < 1 || isempty(n), n = 10; end;

coeffs = wiresaw1(n,v);
K = coeffs{1};
D = coeffs{2};
M = coeffs{3};

K = K + eta*D;
D = D + eta*eye(n);


coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end