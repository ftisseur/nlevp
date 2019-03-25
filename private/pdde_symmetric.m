function [coeffs,fun,F] = pdde_symmetric(n)
%PDDE_SYMMETRIC  n-by-n NEP from a partial delay differential equation.
%  [COEFFS,FUN,F] = nlevp('pdde_symmetric',n) constructs a
%  (n-1)^2-by-(n-1)^2 nonlinear eigenvalue problem 
%  M+A - lambda*I + e^(-2*lambda)*B arising in the study of a partial delay
%  differential equation. The eigenvalues of interest are the ones near zero. 
%  The matrices are returned in a cell array: COEFFS = {M+A, -I, B}.
%  FUN is a function handle to evaluate the functions 1, lambda,
%  e^(-2*lambda) and their derivatives.
%  F is the function handle M+A - lambda*I + e^(-2*lambda)*B.
%  This problem has the properties nep, real, symmetric, scalable, sparse,
%  banded.

%  Reference:
%  MATLAB file written by Fei Xue. March 2019.

if nargin < 1 || isempty(n)
    n = 128;
end
[MSX,MSY] = meshgrid(0:pi/n:pi,0:pi/n:pi);

MSX = MSX(2:end-1,2:end-1); 
MSY = MSY(2:end-1,2:end-1);

da = (sin(MSX).*sin(MSY)).^2; db = (sin(MSX+MSY)+1.31);

A = spdiags(reshape(da,(n-1)^2,1),0,(n-1)^2,(n-1)^2);
B = spdiags(reshape(db,(n-1)^2,1),0,(n-1)^2,(n-1)^2);
M = -delsq(numgrid('S',n+1))/((pi/n)^2);

coeffs{1} = M+A; 
coeffs{2} = -speye(size(coeffs{1})); 
coeffs{3} = B;

fun = @(lam) pddesymm_fun(lam);
F = @(lam) coeffs{1} + lam*coeffs{2} + exp(-2*lam)*coeffs{3};

end


function varargout = pddesymm_fun(lam)
lam = lam(:);
n = length(lam);

varargout{1} = [ones(n,1), lam, exp(-2*lam)];
if nargout >= 2
    f = -2*exp(-2*lam);
    varargout{2} = [zeros(n,1), ones(n,1), f];
    for i = 2 : nargout-1
        f = -2*f;
        varargout{i+1} = [zeros(n,2), f];
    end
end

end