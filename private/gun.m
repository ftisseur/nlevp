function [coeffs,fun,F,xcoeffs] = gun
%GUN  NEP from model of a radio-frequency gun cavity.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('gun') returns 9956-by-9956 matrices
%  K, M, W1 and W2 for the nonlinear eigenvalue problem defined by
%  T(lam)*x = [K-lam*M+i*lam^(1/2)*W1+i*(lam-108.8774^2)^(1/2)*W2]*x = 0.
%  The matrices are returned in a cell array: COEFFS = {K,M,W1,W2}.
%  FUN is a function handle to evaluate the functions 1,-lam,i*lam^(1/2),
%  and i*(lam-108.8774^2)^(1/2) and their derivatives.
%  F is the function handle T(lam).
%  XCOEFFS is the cell {1, 1, L1, L2; K, M, R1', R2'} to exploit the low
%  rank of W1 = L1*R1' and W2 = L2*R2'.
%  This problem has the properties nep, sparse, banded (843 bands),
%  low-rank.

%  References:
%  Ben-Shan Liao, Subspace projection methods for model order reduction and
%  nonlinear eigenvalue computation, PhD Thesis, Department of Mathematics,
%  University of California at Davis, 2007.

load gun

coeffs = {K, M, W1, W2};
fun = @gun_fun;
F = @(lam) K - lam*M + 1i*lam^(.5)*W1+1i*(lam-108.8774^2)^(.5)*W2;
xcoeffs = {1, 1, L1, L2; K, M, R1', R2'};
end

function varargout = gun_fun(lam)

lam = lam(:);
n = length(lam);
sigma2 = 108.8774^2;
f3 = 1i*sqrt(lam);
f4 = 1i*sqrt(lam-sigma2);

varargout{1} = [ones(n,1),-lam,f3,f4];
if nargout >= 2,
    f3 = 0.5*f3./lam;
    f4 = 0.5*f4./(lam-sigma2);
    varargout{2} = [zeros(n,1),-ones(n,1),f3,f4];
end
for i=2:nargout-1,
    f3 = (3/2-i)*f3./lam;
    f4 = (3/2-i)*f4./(lam-sigma2);
    varargout{i+1} = [zeros(n,2),f3,f4];
end


end
