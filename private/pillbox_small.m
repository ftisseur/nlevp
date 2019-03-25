function [coeffs,fun,F] = pillbox_small(problem)
%PILLBOX_SMALL  20-by-20 NEP from a RF pillbox cavity.
%  [COEFFS,FUN,F] = nlevp('pillbox_small','problem') returns a 20-by-20 
%  matrices K, M, W1,and W2 for the nonlinear eigenvalue problem defined by
%  T(lam)*x = [K - lam*M +i*(lam-sigma1^2)^(1/2)*W1 + 
%  i*(lam-sigma^2)^(1/2)*W3]*x = 0.
%  The variable PROBLEM can take the value 'a' (default), 'b' or 'brd',
%  and allows variations of the same problem.
%  COEFFS returns the matrices {K,M,W1,W2} in a cell array.
%  FUN is a function handle to evaluate the functions 1, -lam,
%  i*(lam-sigma1^2)^(1/2) and i*(lam-sigma2^2)^(1/2) and their derivatives.
%  F is the function handle T(lam).
%  The problem has the property nep, symmetric.

%  Reference:
%  https://github.com/oamarques/PAL

%  R. Van Beeumen, Osni Marques et al, Computing resonant modes
%  of accelerator cavities by solving nonlinear eigenvalue problems
%  via rational approximation, J. Comput. Phys., 374:1031-1043, 2018.

if nargin < 1
    problem = 'a';
end

if problem == 'a'
    load('pillbox_smallA.mat')
    sigma1 =  1.5000e-01;
    sigma2 = 0.7;
elseif problem == 'b'
    load('pillbox_smallB.mat')
    sigma1 =  1.2;
    sigma2 = 2.4;
elseif strcmp(problem, 'brd')
    load('pillbox_smallBRD.mat')
    sigma1 =  1.2;
    sigma2 = 2.4;
else % wrong input
    fprintf('This problem does not exist. Call  nlevp(help, ''pillbox_small'') for more information.')
    error
end

coeffs = {K, M, W1, W2};
fun = @(lam) pillbox_fun(lam, sigma1, sigma2);
F = @(lam) K -lam*M + 1i*(lam -sigma1^2)^0.5*W1 + 1i*(lam-sigma2^2)^0.5*W2;

end

function varargout = pillbox_fun(lam, sigma1, sigma2)

lam = lam(:);
n = length(lam);

% cutoff frequencies


f1 = 1i*sqrt(lam - sigma1^2);
f2 = 1i*sqrt(lam - sigma2^2);

varargout{1} = [ones(n,1),-lam,f1,f2];
if nargout >= 2
    
    f1 = 0.5*f1./(lam-sigma1^2);
    f2 = 0.5*f2./(lam-sigma2^2);
    varargout{2} = [zeros(n,1),-ones(n,1),f1,f2];
end
for i=2:nargout-1
    f1 = (3/2-i)*f1./(lam-k1);
    f2 = (3/2-i)*f2./(lam-k3);
    varargout{i+1} = [zeros(n,2),f1,f2];
end

end
