function [coeffs,fun,F] = surveillance
%SURVEILLANCE  21-by-16 QEP from surveillance camera callibration.
%  [COEFFS,FUN,F] = nlevp('surveillance') generates a 21-by-16 quadratic
%  matrix polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0
%  arising from calibration of a surveillance camera using a human body as a
%  calibration target.  The eigenvalue represents the focal length of the
%  camera.  This particular data set is synthetic and corresponds to a
%  600-by-400 pixel camera
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2,
%  and their derivatives.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  This problem has the properties pep, qep, real, nonsquare.

%  Reference:
%  B. Micusik and T. Pajdla. Simultaneous surveillance camera
%  calibration and foot-head homology estimation from human
%  detections. In IEEE Computer Society Conference on Computer Vision
%  and Pattern Recognition (CVPR), San Francisco, USA, 2010.

load surveillance

% OUTPUT
coeffs = {A0, A1, A2};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
