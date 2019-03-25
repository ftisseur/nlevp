function [coeffs,fun,F] = wiresaw1(n,v)
%WIRESAW1   Gyroscopic QEP from vibration analysis of a wiresaw.
%  [COEFFS,FUN,F] = nlevp('wiresaw1',N,V) constructs the N-by-N gyroscopic
%  quadratic matrix polynomial lambda^2*M + lambda*D + K arising
%  from the vibration analysis of a wiresaw.
%  V is a nonnegative parameter corresponding to the speed of the wire.
%  The default values are N = 10 and V = 0.01.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2,
%  and their derivatives.
%  F is the function handle lambda^2*M + lambda*D + K.
%  This problem has the properties pep, qep, t-even, gyroscopic,
%  parameter_dependent, scalable.

%  Reference:
%  S. Wei and I. Kao, Vibration analysis of wire and frequency response in
%  the modern wiresaw manufacturing process, Journal of Sound and Vibration
%  213(5) (2000), pp. 1383-1395.

if nargin < 2 || isempty(v), v = 0.01; end;
if nargin < 1 || isempty(n), n = 10; end;

M = eye(n);
D = zeros(n);
K = zeros(n);
for j = 1:n
    for k = 1:n
        if rem(j+k,2), D(j,k) = 8*j*k*v/(j^2-k^2); end
    end
    K(j,j) = j^2*pi^2*(1-v^2);
end


coeffs = {K/2, D/2, M/2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end