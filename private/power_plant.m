function [coeffs,fun,F] = power_plant(mu)
%POWER_PLANT   8-by-8 QEP from simplified nuclear power plant problem.
%[COEFFS,FUN,F] = nlevp('power_plant',MU) generates the matrices of an 8-by-8
%  quadratic matrix polynomial lambda^2*M+lambda*D+K arising in the study of
%  the dynamic behaviour of a nuclear power plant. The matrices D and M are
%  real symmetric and K is defined as K = (1+i*MU)K0 with K0 real symmetric.
%  The parameter MU describes hysteretic damping added to the model.
%  The default value for MU is 0.2.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*M+lambda*D+K.
%  This problem has the properties pep, qep, symmetric, parameter_dependent.

%  Reference:
%  T. Itoh, Damped vibration mode superposition method for dynamic response
%     analysis, Earthquake Eng. Struct. Dyn. 2 (1973), pp. 45--57
%  F. Tisseur and K. Meerbergen, The quadratic eigenvalue problem, SIAM
%     Rev. 43 (2001), pp. 235--286

if nargin < 1 || isempty(mu), mu = 0.2; end

M = diag([0.54e1 0.168e7 0.552e2 0.8711e8 0.578e2 0.235e9 0.16e2 0.214e8]);

D = [0  -0.7492e4 -0.999e1 0.979e4 0 0 0 0
     0 0 0.7492e4 -0.1221e8 0 0 0 0
     0 0 0 -0.9674e6 -0.538e3 -0.538e6 0 0
     0 0 0 0 0.9576e6 -0.1652e10 0 0
     0 0 0 0 0 -0.5042e7 -0.3e2 0.57e5
     0 0 0 0 0 0  -0.3e5 -0.182e9
     0 0 0 0 0 0 0 -0.57e5
     0 0 0 0 0 0 0 0];

DD = diag([0.999e1 0.1045e8 0.548e3 0.4329e10 0.6178e4 ...
           0.4343e11 0.3e2 0.3473e9]);
D = D+D'+DD;

K = [0 -0.87e7 -0.116e5 0.1137e8 0 0 0 0
     0 0 0.87e7 -0.1723e11 0 0 0 0
     0 0 0 -0.9334e9 -0.518e6 -0.518e9 0 0
     0 0 0 0 0.922e9 -0.6878e13 0 0
     0 0 0 0 0 0.3633e9 -0.617e4 0.1172e8
     0 0 0 0 0 0 -0.617e7 -0.2553e12
     0 0 0 0 0 0 0 -0.1172e8
     0 0 0 0 0 0 0 0];

DD = diag([0.116e5 0.1522e11 0.5296e6 0.9461e13 0.6862e6 0.9953e13 ...
           0.6170e4 0.2893e12]);
K = K+K' + DD;
K = (1+sqrt(-1)*mu)*K;

coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
