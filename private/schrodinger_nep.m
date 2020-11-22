function [coeffs, fun, F] = schrodinger_nep(n, alpha)
%SCHRODINGERnep   N-by-N NEP from a Schroedinger equation
% [coeffs, fun, F] = nlevp('schrodingernep', N, alpha) generates the N-by-N
% nonlinear matrix function A - lam*I + g(lam)*B1 + f(lam)*B2, which
% iderives from the discretization of a boundary value problem from a
% Schroedinger equation. The potential function V(z) is defined as
% 1+sin(alpha*z) in [0,1] and is a constant equal to 10 in ]1,8[. The
% scalar functions are g(lam) = cosh(7sqrt(lam+V0)), and 
% f(lam) = sinh(7sqrt(lam+V0))/sqrt(z+V0).
% The matrices are defined as follows: A = D + W, where D is the classic
% second difference operator on n points, but with D(n,n-1:n) = [0 0], and
% W is diag(V(1/n), ..., V(1-1/n), 0); I is the identity matrix, but 
% I(n,n) = 0; B1 is e_n*e_n^T and B2 = n/2*e_n*e_{n-2}^T - 2n*e_n*e_{n-1}^T
% + 3n/2*e_n*e_n^T.
% COEFFS is the cell of the matrix coefficients {A, I, B1, B2},
% FUN is a function handle to evaluate the functions [1, -lam,
% g(lam), f(lam)].
% F is a function handle to evaluate the matrix-valued function A - lam*I +
% g(lam)*B1 + f(lam)*B2.
% The problem has the properties nep, real, scalable, sparse,
% parameter-dependent, tridiagonal, banded.

% References: Problem from  NEP-Pack Library: "NEP-PACK: A Julia package
% for nonlinear eigenproblems", E. Jarlebring and M. Bennedich and G. Mele 
% and E. Ringh and P. Upadhyaya, 2018, arXiv:1811.09592. Available at
% https://github.com/nep-pack.

if nargin < 1 || isempty(n)
    n = 20;
end
if nargin < n || isempty(alpha)
    alpha  = 25*pi/2;
end

h = 1/n;
I = speye(n);
I(n,n) = 0;
e = ones(n,1);
D = spdiags([e -2*e e], -1:1, n, n);
D(n,n-1) = 0; D(n,n) = 0;
D = D*n^2;
Vvec = [1 + sin(alpha*[1:n+1]'*h); 0];
A = D - spdiags(Vvec, 0, n,n);
F = sparse(n,n);
F(end, end-2) = n/2;
F(end, end-1) = -2*n;
F(end, end) = 3*n/2;

% Setting the coefficients
coeffs{1} = A;
coeffs{2} = I;
coeffs{3} = sparse(n,n);
coeffs{3}(end,end) = 1;
coeffs{4} = F;

% setting the functions
L0 = 1;
L1 = 8;
V0 = 10;
f = @(lam) sinh((L1-L0).*sqrt(lam+V0))./sqrt(lam+V0);
g = @(lam) cosh((L1-L0).*sqrt(lam+V0));
F = @(lam) coeffs{1} - lam*coeffs{2} + g(lam)*coeffs{3} +f(lam)*coeffs{4};
fun = @(lam) schrodinger_nep_fun(lam, L1, L0, V0);

end


function varargout = schrodinger_nep_fun(lam, L1, L0, V0)
% We do not return any derivative because they are too complicated to
% compute

lam = lam(:);
l = length(lam);
f = @(lam) sinh((L1-L0).*sqrt(lam+V0))./sqrt(lam+V0);
g = @(lam) cosh((L1-L0).*sqrt(lam+V0));
varargout{1} = [ones(l,1) -lam, g(lam), f(lam) ];

end