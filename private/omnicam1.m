function [coeffs,fun,F] = omnicam1
%OMNICAM1  9-by-9 QEP from model of omnidirectional camera.
%  [COEFFS,FUN,F] = nlevp('omnicam1') generates a 9-by-9 quadratic matrix 
%  polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0 arising from a
%  model of an omnidirectional camera (one with angle of view greater
%  than 180 degrees).
%  A0 has rank 1, A2 has rank 5, and A1 has full rank,
%  in each case the rank being the number of nonzero columns.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*A2 + lambda*A1 + A0.
%  This problem has the properties pep, qep, real.

%  Reference:
%  B. Micusik and T. Pajdla. Estimation of omnidirectional camera
%  model from epipolar geometry. In  IEEE Computer Society Conference
%  on Computer Vision and Pattern Recognition (CVPR'03). IEEE Computer
%  Society, Los Alamitos, CA, USA, 2003.

load omnicam1

coeffs = {D3, D2, D1};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
