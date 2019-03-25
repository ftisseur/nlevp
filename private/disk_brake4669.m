function [coeffs,fun,F] = disk_brake4669(omega)
%DISK_BRAKE4669  4669-by-4669 QEP from a disk brake model.
%  [COEFFS,FUN,F] = nlevp('disk_break',omega) constructs a 4669-by-4669
%  quadratic matrix polynomial lambda^2*M + lambda*D(omega) + K(omega)
%  arising in the study of a disk brake model. The matrices are
%  D(omega) = D1 + DR/omega + omega*DG + D4/Fref,
%  K(omega) = K1 + KR + omega^2*KGEO,
%  where Fref = 1600 and omega is a real scalar (e.g., 2pi <= omega <= 8pi).
%  The matrices are returned in a cell array: COEFFS = {K(omega), D(omega), M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K(omega) + lambda*D(omega) + lambda^2*M.
%  This problem has the properties pep, qep, real, parameter_dependent,
%  sparse

%  Reference:
%  N. Graebner, V. Mehrmann, S. Quraishi, C. Schroeder and U. von Wagner.
%  Numerical methods for parametric model reduction in the simulation of
%  disc brake squeal. Z. Angew. Math. Mech., 96(12):1388-1405, 2016.

if nargin < 1 || isempty(omega)
    omega = 2*pi;
end
if ~isreal(omega), error('omega must be real');end

Fref = 1600;
omegaref = 1;
load disk_brake4669;
D = D1 + DR/omega + omega*DG + D4/Fref;
K = K1 + KR + (omega/omegaref)^2*KGEO;
coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
