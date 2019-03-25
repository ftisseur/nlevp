function [coeffs,fun,F,xcoeffs] = acoustic_wave_2d(n,z)
%ACOUSTIC_WAVE_2D   Acoustic wave problem in 2 dimensions.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('acoustic_wave_2d',N,Z) constructs a
%  quadratic matrix polynomial lambda^2*M + lambda*D + K that arises from
%  the discretization by finite elements of a two-dimensional acoustic
%  wave equation.  The coefficient matrices are M-by-M with
%  M = N1*(N1-1), H = 1/N1 being the mesh size, with N1 chosen such
%  that M approximates N.  The damping matrix has the
%  form 2*pi*i*Z^(-1)*C, where C = H*L*L' is a low rank real symmetric
%  matrix and the scalar parameter Z is the impedance (possibly complex).
%  The default values are N = 30 (N1 = 6) and Z = 1.
%  The eigenvalues lie in the upper half of the complex plane.
%  The matrices are returned in a cell array: COEFFS = {K, D, M};
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  XCOEFFS returns the cell {1 2*pi*i*L 1;K L' M} to exploit the low rank
%  of D. 
%  This problem has the properties pep, qep, symmetric, *-even
%  parameter-dependent, scalable, sparse, banded (6 bands), low-rank.

%  Reference:
%  F. Chaitin-Chatelin and M. B. van Gijzen, Analysis of parameterized
%  quadratic eigenvalue problems in computational acoustics with homotopic
%  deviation theory, Numer. Linear Algebra Appl. 13 (2006), pp. 487-512

if nargin < 2 || isempty(z)
    z = 1; 
end

if nargin < 1 || isempty(n)
    n = 30;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
            'now targets the true dimension of the problem'])
end

% Find out best n1 with n1*(n1-1) approx n.
n1 = floor(0.5 + sqrt(n+0.25));
if abs(n-n1*(n1-1)) > abs(n-(n1+1)*n1)
    n1 = n1 + 1;
end
n1 = max(n1,2);
% n = n1*(n1-1)

h = 1/n1;
D = spdiags(ones(n1,1)*[-1,4,-1],-1:1,n1,n1);
D(n1,n1) = 2;
T = spdiags(ones(n1-1,2),[-1,1],n1-1,n1-1);
S = speye(n1); 
S(n1,n1) = 0.5;
E = sparse(n1,n1,1,n1,n1);

M = h^2*kron(speye(n1-1),S);
C = (h/z)*kron(speye(n1-1),E);
K = kron(speye(n1-1),D) + kron(T,-S);

coeffs = {K, 2*pi*1i*C, -(2*pi)^2*M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);
% Building xcoeffs
en = sparse(zeros(n1,1));
en(n1) = 1;
L = kron(eye(n1-1),en);
xcoeffs1 = {1, 2*pi*1i*h*L/z, 1};
xcoeffs2 = {coeffs{1}, L', coeffs{3}};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};

end
