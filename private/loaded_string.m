function [coeffs,fun,F,xcoeffs] = loaded_string(n,kappa,m)
%LOADED_STRING    REP from finite element model of a loaded vibrating string.
%   [COEFFS,FUN,F,XCOEFFS] = nlevp('loaded_string',N,KAPPA,m) generates an
%   N-by-N rational matrix function A - lambda*B + lambda/(lambda-sigma)*C
%   with sigma = KAPPA/m arising in the finite element discretization
%   of a boundary problem describing the eigenvibration of a string
%   with a load of mass m attached by an elastic spring of stiffness KAPPA.
%   The default values are N = 20, KAPPA = 1, m=1.
%   The matrices are returned in a cell array: COEFFS = {A, B, C}.
%   FUN is a function handle to evaluate the functions 1,-lambda, and
%   lambda/(lambda-sigma) and their derivatives.
%   F is the function handle A - lambda*B + lambda/(lambda-sigma)*C.
%   XCOEFFS is the cell {1, 1, en; A, B, enz'}, where C = en*enz', 
%   that exploits the 1-rank structure of C. 
%   This problem has the properties rep, qep, real, symmetric,
%   parameter_dependent, scalable, sparse, tridiagonal, sparse, low-rank.

%   Reference:
%   S. I. Solov"ev. Preconditioned iterative methods for a class of nonlinear
%   eigenvalue problems. Linear Algebra Appl., 415 (2006), pp. 210-229.

if nargin < 1 || isempty(n), n = 20; end;
if nargin < 2 || isempty(kappa), kappa = 1; end;
if nargin < 3 || isempty(m), m = 1; end;

% A = toeplitz([2,-1,zeros(1,n-2)]); A(n,n) = 1;
A = spdiags(ones(n,1)*[-1,2,-1],-1:1,n,n); A(n,n) = 1;
% B = toeplitz([4,1,zeros(1,n-2)]); B(n,n) = 2;
B = spdiags(ones(n,1)*[1,4,1],-1:1,n,n); B(n,n) = 2;
% C = zeros(n); C(n,n) = 1;
C = sparse(n,n,1,n,n);

coeffs = {n*A, 1/(6*n)*B, kappa*C};
fun = @(lam) loaded_string_fun(lam,kappa/m);
F = @(lam) coeffs{1} - lam*coeffs{2} + lam/(lam-kappa/m)*coeffs{3}; 
%Building xcoeffs
en = sparse(zeros(n,1)); 
en(n) = 1;
enz = en;
enz(n) = kappa;
xcoeffs = {1, 1, enz; coeffs{1}, coeffs{2}, en'};
end


function varargout = loaded_string_fun(lam,sigma)

lam = lam(:);
n = length(lam);

varargout{1} = [ones(n,1),-lam,lam./(lam-sigma)];
if nargout >= 2
    f = -sigma./((lam-sigma).^2);
    varargout{2} = [zeros(n,1),-ones(n,1),f];
    for i=2:nargout-1
        f = -i*f./(lam-sigma);
        varargout{i+1} = [zeros(n,2),f];
    end
end
end
