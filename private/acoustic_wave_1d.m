function [coeffs,fun,F,xcoeffs] = acoustic_wave_1d(n,z)
%ACOUSTIC_WAVE_1D   Acoustic wave problem in 1 dimension.
%  [COEFFS,FUN,F] = nlevp('acoustic_wave_1d',N,Z) constructs an N-by-N
%  quadratic matrix polynomial lambda^2*M + lambda*D + K that arises
%  from the discretization of a 1D acoustic wave equation.
%  The damping matrix has the form 2*pi*i*Z^(-1)*C, where
%  C = e_n*e_n', where e_n = [0 ... 0 1]', and the scalar parameter Z is
%  the impedance (possibly complex).
%  The default values are N = 10 and Z = 1.
%  The eigenvalues lie in the upper half of the complex plane.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  XCOEFFS returns the cell {1 en 1;K en' M} to exploit the low rank of D.
%  This problem has the properties pep, qep, symmetric, *-even
%  parameter-dependent, scalable, sparse, tridiagonal, banded, low-rank.

%  Reference:
%  F. Chaitin-Chatelin and M. B. van Gijzen, Analysis of parameterized
%  quadratic eigenvalue problems in computational acoustics with homotopic
%  deviation theory, Numer. Linear Algebra Appl. 13 (2006), pp. 487-512.

if nargin < 2 || isempty(z)
    z = 1; 
end
if nargin < 1 || isempty(n) 
    n = 10; 
end

h = 1/n;

e = ones(n,1);
K = spdiags([-e,2*e,-e],-1:1,n,n);
K(n,n) = 1;
K = n*K;

D = sparse(n,n,1/z,n,n);
M = speye(n); M(n,n) = 0.5; M = h*M;

coeffs = {K, 2*pi*1i*D, -(2*pi)^2*M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);
en = sparse(zeros(n,1)); en(n) = 1;
enz = en; enz(n) = 2*pi*1i/z;
xcoeffs1 = {1,  enz, 1};
xcoeffs2 = {coeffs{1}, en', coeffs{3}};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};
end
