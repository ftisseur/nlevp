function [coeffs,fun,F] = damped_gyro(n,C,epsi)
%DAMPED_GYRO  QEP from a damped gyroscopic system.
%  [COEFFS,FUN,F] = nlevp('damped_gyro',N,C,epsi) builds a quadratic
%  matrix polynomial lambda^2*M + lambda*(G + epsi*D) + K that arises from
%  the discretization by finite differences of a damped gyroscopic
%  system. The coefficient matrices M = M^T > 0, K = K^T >0, D = D^T > 0,
%  G = -G^T are N1-by-N1, where N1 = (ceil(sqrt(N)))^2. Every
%  matrix M, G, K is the sum of two scaled kronecker products, for
%  example M = kron(c_11*I, M1) + kron(c_12*M1, I). The parameter C is
%  a 4-by-2 matrix that contains the scalars c_ij > 0. The default values 
%  are N = 36, epsi = 10e-4 and C = [1 1.3; 0.1 1.1; 1 1.2;1.05 0.9]. 
%  If epsi = 0, only the first three rows of C are used.
%  The matrices are returned in a cell array: COEFFS = {K, G+epsi*D, M};
%  FUN is a function handle to evaluate the monomials
%  1,lambda,lambda^2 and their derivatives.
%  F is the function handle K + lambda*(G+epsi*D) + lambda^2*M.
%  This problem has the properties pep, qep, real, parameter-dependent,
%  scalable, sparse, banded (7 bands).

%  Reference: Problem 5.3 (and 5.1) from
%  T.-M. Hwang, W.-W. Lin, V. and Mehrmann. Numerical solution of quadratic
%  eigenvalue problems with structure-preserving methods.
%  SIAM J. Sci. Comput., 24(4):1283-1302, 2003.

if nargin < 1 || isempty(n)
    n = 36;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
            ' targets the true dimension of the problem'])
end;
if nargin < 2 || isempty(C)
   C = [1 1.3; 0.1 1.2; 1 1.2; 1.05 0.9];
end;
if nargin < 3 || isempty(epsi)
   epsi = 1e-4;
end
if ~isequal([4,2], size(C))
    error('C must be either empty or a 4-by-2 matrix')
end

m = ceil(n^0.5);
n1 = (m)^2;
n1 = max(n1,2);
I = speye(m);
B1 = spdiags(ones(m-1,1),-1,m,m);
G1 = B1 - B1';
M1 = (4*I + B1 + B1')/6;
K1 = B1 + B1' - 2*I;
D1 = B1 + B1' + 2*I;

M = kron(C(1,1)*I, M1) + kron(C(1,2)*M1, I);
G = kron(C(2,1)*I, G1) + kron(C(2,2)*G1, I);
K = kron(C(3,1)*I, K1) + kron(C(3,2)*K1, I);

if epsi ~= 0
    D = kron(C(4,1)*I, D1) + kron(C(4,2)*D1, I);
    coeffs = {K, G+epsi*D, M};
else
    coeffs = {K, G, M};
end

fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
