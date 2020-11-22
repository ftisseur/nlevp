function [coeffs,fun,F,xcoeffs] = schrodinger_abc(n)
%SCHRODINGER_ABC  NEP from Schrodinger equation with absorbing boundary condition.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('schrodinger_abc',N) returns N-by-N (N>2)
%  matrices A, B, G and F for the nonlinear eigenvalue problem defined by
%  [A-lambda*B+g(lambda)*G+h(lambda)*H]*x = 0, where
%  g(lambda) = cosh((L1-L0)*sqrt(lambda+V0)) and
%  f(lambda) = sinh((L1-L0)*sqrt(lambda+V0))/sqrt(lambda+V0)
%  with L0, L1, V0 real constants.
%  The matrices are returned in a cell array: COEFFS = {A, B, G, H}.
%  FUN is a function handle to evaluate the functions 1, -lambda, g(lambda),
%  and h(lambda).
%  F is the function handle F(lambda) = A-lambda*B+g(lambda)*G+h(lambda)*H.
%  XCOEFFS is the cell {1, 1, en, en; A, B, en', u} to exploit the low rank
%  of G and H.
%
%  This problem has the properties nep, real, sparse, banded,
%  low-rank, scalable.

% References: Problem taken from Tutorial 1 (ABC), NEP-Pack Library:
% "NEP-PACK: A Julia package for nonlinear eigenproblems", E. Jarlebring
% and M. Bennedich and G. Mele and E. Ringh and P. Upadhyaya, 2018,
% arXiv:1811.09592. Available at https://github.com/nep-pack.


if nargin < 1 || isempty(n), n = 1000; end
if n<3, error('n must be larger than 2');end

L0 = 1; L1 = 8; V0 = 10;
h = 1/(n-1);
xv = (0:h:L0);
alpha = 25*pi/2;
V = 1+sin(alpha*xv); V(end) = 0;
hh = 1/h^2;

A = diag(sparse(hh*[ones(n-2,1);0]),-1)+ diag(sparse(hh*[-2*ones(n-1,1);0]))+...
    diag(sparse(hh*ones(n-1,1)),1);
A = A-diag(sparse(V));
B = diag(sparse([ones(n-1,1);0]));
G = sparse([n],[n],[1]);
H = sparse([n, n, n],[n-2, n-1, n],[1/(2*h), -2/h, 3/(2*h)]);

coeffs = {A,B,G,H};

f = @(lam) sqrt(lam+V0);
g1 = @(lam) cosh((L1-L0)*f(lam));
g2 = @(lam) sinh((L1-L0)*f(lam))./f(lam);

fun = @(lam) schrodinger_abc_fun(lam,L0,L1,V0);
F = @(lam) coeffs{1}-lam*coeffs{2}+coeffs{3}*g1(lam)+coeffs{4}*g2(lam);

en = sparse(zeros(n,1)); en(n) = 1;
enaux = en;
enaux(end-2:end) = [1/(2*h); -2/h; 3/(2*h)];

xcoeffs1 = {1,  1, en, en};
xcoeffs2 = {coeffs{1}, coeffs{2}, en', enaux'};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};

end


function varargout = schrodinger_abc_fun(lam,L0,L1,V0)
% We do not return any derivatives because they are too complicated to
% compute

lam = lam(:);
f = @(lam) sqrt(lam+V0);
g1 = @(lam) cosh((L1-L0)*f(lam));
g2 = @(lam) sinh((L1-L0)*f(lam))./f(lam);
%First derivatives
%gp = .5*(L1-L0)*h(lam);
%hp = .5*(L1-L0)*g(lam)./(lam+V0) - .5*h(lam)./(lam+V0);

n = length(lam);
varargout{1} = [ones(n,1),-lam, g1(lam), g2(lam)];
%fp = [zeros(n,1),-ones(n,1),gp(lam), fp(lam)];
end



