function [coeffs,fun,F] = disk_brake100(omega)
%DISK_BRAKE100 100-by-100 QEP from a disk brake model.
%  [COEFFS,FUN,F] = nlevp('disk_break',omega) constructs a 100-by-100
%  quadratic matrix polynomial lambda^2*M + lambda*D(omega) + K
%  arising in the study of a disk brake model. The matrices are
%  D(omega) = D1 + DR/omega and K = K1 + KR, with omega a real scalar 
%  (e.g., 2pi <= omega <= 8pi).
%  The matrices are returned in a cell array: COEFFS = {K, D(omega),M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D(omega) + lambda^2*M.
%  This problem has the properties pep, qep, real, parameter_dependent.

%  Reference:
%  N. Graebner, V. Mehrmann, S. Quraishi, C. Schroeder and U. von Wagner.
%  Numerical methods for parametric model reduction in the simulation of 
%  disc brake squeal. Z. Angew. Math. Mech., 96(12):1388-1405, 2016.

if nargin < 1 || isempty(omega)
    omega = 2*pi; 
end

if~isreal(omega), error('omega must be real');end 

load disk_brake100;

D = D1 + DR/omega + omega*DG;
K = K1 + KR;
coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end