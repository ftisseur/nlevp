function [coeffs,fun,F] = sign1(n,a)
%SIGN1  QEP from rank-1 perturbation of sign operator.
%  [COEFFS,FUN,F] = nlevp('sign1',N,MU) generates the coefficient matrices of
%  a quadratic matrix polynomial lambda^2*A + lambda*B + C of dimension N.
%  N should be odd. Otherwise N+1 is used.
%  The spectrum of this matrix polynomial is the second order
%  spectrum of a rank one perturbation of the operator of multiplication
%  by the function sign(x) acting in L^2(-\pi,\pi), with respect to
%  the Fourier basis: e^{-iNx},\ldots,1,\ldots e^{iNx}.
%  The strength of the perturbation is the real parameter MU.
%  By default N = 81 and MU = 0.
%  The matrices are returned in a cell array: COEFFS = {C, B, A}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*A + lambda*B + C.
%  This problem has the properties qep, hermitian, parameter-dependent,
%  scalable.

%  Reference:
%  L. Boulton, Non-variational approximation of discrete eigenvalues of
%  self-adjoint operators, IMA J. Numer. Anal. 27 (2007), pp. 102-121.

if nargin < 1 || isempty(n)
    n = 81;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
            'now targets the true dimension of the problem']);
end
if nargin < 2, a = 0; end

m = round((n-1)/2);
n = 2*m+1;

% Fourier coeff. of the symbol.
fc = zeros(n,2);
fc(1,2) = 1;                        % FC of the square of the operator
for k=2:2:2*m
    fc(k,1) = 1i*2/(pi*(k-1));      % FC of the operator
end

% Matrix truncations
K = zeros(n);
K(m+1,m+1) = a;
MK = zeros(n);
B = toeplitz(fc(1:n,1));
MK(:,m+1) = B(:,m+1);
Q = toeplitz(fc(1:n,2));
A0 = Q + a*MK + a*MK' + K^2;
A1 = -2*(B+K);
A2 = eye(n);

% OUTPUT
coeffs = {A0, A1, A2};
fun = @(lam)nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
