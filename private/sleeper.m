function [coeffs,fun,F,sol] = sleeper(n)
%SLEEPER   QEP modelling a railtrack resting on sleepers.
%  [COEFFS,FUN,F] = nlevp('sleeper',N) or
%  [COEFFS,FUN,F,SOL] = nlevp('sleeper',N), where N >= 5, generates
%  an N-by-N proportionally damped quadratic matrix polynomial
%  Q(lambda) = lambda^2 A2 + lambda*A1 + A0
%  arising from the analysis of a railtrack resting on sleepers.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2 A2 + lambda*A1 + A0.
%  The exact eigenvalues and eigenvectors are known and are returned as
%  SOL.EVAL and SOL.EVEC. Note that A0,A1,A2 are sparse, but SOL.EVEC is
%  full. SOL.EVEC is constructed only if SOL is requested.
%  The default dimension of the problem is N = 10.
%  The problem has  properties pep, qep, real, scalable,
%  proportionally-damped, symmetric, sparse, solution.

%  Reference:
%  P. Lancaster, and P. Rozsa, The spectrum and stability of a vibrating
%  rail supported by sleepers, Computers Math. Applic., 31(4/5):201-213,
%  1996

if nargin < 1, n = 10; end

if n < 5, error('N must be at least 5.'), end
% Solution constructed below is incorrect for N < 5.

% Create circulant matrix A
% A = toeplitz([-2 1 zeros(1,n-2)]);
A = spdiags(ones(n,1)*[1,-2,1],-1:1,n,n);
A(1,n) = 1; A(n,1) = 1;

% AA is the circulant matrix A^2.
% AA = toeplitz([6 -4 1 zeros(1,n-3)]);
AA = spdiags(ones(n,1)*[1,-4,6,-4,1],-2:2,n,n);
AA(1,n) = -4; AA(n,1) = -4; AA(1,n-1) = 1; AA(2,n) = 1;
AA(n-1,1) = 1; AA(n,2) = 1;

I = speye(n);

%OUTPUT
coeffs = {I + A + AA, I + AA, I};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

% Generate exact eigenvalues and eigenvectors, if asked.
if nargout > 3
    sol.evec = 1/sqrt(n)*fft(eye(n));
    sol.evec = [sol.evec,sol.evec];
    mu = -4*sin(pi*(0:n-1)'/n).^2;
    l1 = -.5*(1+mu.^2+sqrt(mu.^4-2*mu.^2-4*mu-3));
    l2 = (1+mu+mu.^2)./l1;
    sol.eval = [l1;l2];
end


end
