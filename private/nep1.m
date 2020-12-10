function [coeffs,fun,F] = nep1
%NEP1  2-by-2 basic NEP example.
%  [COEFFS,FUN,F] = nlevp('nep1') generates a 2-by-2
%  nonlinear matrix function A + e^(i*lambda^2)*B.
%  This problem has a double non-semisimple eigenvalue in lambda = 0. 
%  All the other eigenvalues, lambda = +-sqrt(2*pi*k) for k integer \neq 0.
%  The matrices are returned in a cell array: COEFFS = {A, B}.
%  FUN is a function handle to evaluate the functions 1, and
%  exp(-i*lambda^2) and their derivatives.
%  F is the function handle  A + e^(i*lambda^2)*B.
%  This problem has the property nep.

%  Reference: 
%  S. Güttel and F. Tisseur. The nonlinear eigenvalue problem. 
%  Acta Numer., 26:1–94, 2017.

n = 2;
A = ones(n,n);
A(1,1) = 0;
B = zeros(n,n);
B(1,1) = 1;

coeffs = {A, B};
fun = @(lam) nep1_fun(lam);
F = @(lam) [ exp(1i*lam^2) 1; 1 1];

end

function varargout = nep1_fun(lam)

lam = lam(:);
n = length(lam);

varargout{1} = [ones(n,1), exp(1i*lam.^2)];
if nargout >= 2
    f = 2*1i*lam.*exp(1i*lam.^2);
    varargout{2} = [zeros(n,1),f];
    for k = 2:nargout-1
        f = 1i.*lam.*f;
        varargout{k+1} = [zeros(n,2),f];
    end
end
end
