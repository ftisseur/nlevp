function [coeffs, fun, F] = square_root()
%SQUARE_ROOT   Square root of a skew-symmetric matrix.
%   [COEFFS,FUN,F] = nlevp('square_root') generates a 20-by-20
%   nonlinear matrix function F(lambda) = A -lambda^(1/2)*I. The matrix A
%   is skew-symmetric and its eigenvalues lie in [4-40i, 4+40i]. The
%   eigenvalues of F(z) are their square roots.
%   This problem has the properties nep, real, sparse. 

%  Reference: S. Elsworth and S. Güttel, "Conversions between barycentric,
%  RKFUN, and Newton representations of rational interpolants", Linear
%  Algebra Appl., 576 (2019), pp. 246–257.


A = 10*gallery('tridiag', 10); 
S = 4*speye(10);
A = [ S , A ; -A , S ];
coeffs = cell(1,2);
coeffs{1} = A;
coeffs{2} = speye(20);
fun = @(lam) square_root_fun(lam);
F = @(lam) coeffs{1} - sqrt(lam)*coeffs{2};

end





function varargout = square_root_fun(lam)
lam = lam(:);
l = length(lam);
varargout{1} = [ones(l,1), -lam.^0.5];
for j = 2:nargout
    varargout{j} = [zeros(l,1), -(0.5-j+2)*lam.^(0.5-j+1)];
end

end
