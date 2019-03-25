function [coeffs,fun,F] = spring(n,mu,tau1,kappa1,tau2,kappa2)
%SPRING    QEP from finite element model of damped mass-spring system.
%  [COEFFS,FUN,F] = nlevp('spring',N,MU,TAU1,KAPPA1,TAU2,KAPPA2) generates
%  an N-by-N quadratic matrix polynomial lambda^2*M + lambda*D + K,
%  where N >= 2, arising from a linearly damped mass-spring system.
%  MU is an N-vector of masses.
%  TAU1 and KAPPA1 are N-vectors of damping constants for the dampers
%  and springs connecting the masses to the ground.
%  TAU2 and KAPPA2 are (N-1)-vectors of damping constants for the dampers
%  and springs connecting adjacent masses.
%  If any of TAU1, TAU2, KAPPA1 or KAPPA2 is a scalar then it is expanded to
%  a vector all of whose elements have that value and TAU1 and KAPPA1
%  then have their first and last elements doubled.
%  The default values are N = 5, MU = 1, TAU1 = 10, KAPPA1 = 5.
%  TAU2 = TAU1(1), KAPPA2 = KAPPA1(1).
%  The matrices are returned in a cell array: coeffs = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*M + lambda*D + K.
%  This problem has the properties pep, qep, real, symmetric,
%  proportionally-damped, parameter_dependent, scalable, sparse,
%  tridiagonal, banded.

%  With default parameters, this function generates the same
%  problem as STRING in the original NELVP release (1.0).

%  Reference:
%  F. Tisseur, Backward error and condition of polynomial
%  eigenvalue problems. Linear Algebra and Appl. 309 (2000), pp. 339-361.

if nargin < 1 || isempty(n), n = 5; end;
if n <= 1, error('N must be at least 2.'), end

if nargin < 2 || isempty(mu), mu = 1; end;
if length(mu) == 1, mu = mu*ones(n,1); end

if nargin < 3 || isempty(tau1), tau1 = 10; end;
if length(tau1) == 1
   tau1 = tau1*ones(n,1); tau1(1) = 2*tau1(1); tau1(n) = 2*tau1(n);
end

if nargin < 4 || isempty(kappa1), kappa1 = 5; end;
if length(kappa1) == 1
   kappa1 = kappa1*ones(n,1); kappa1(1) = 2*kappa1(1); kappa1(n) = 2*kappa1(n);
end

if nargin < 5 || isempty(tau2), tau2 = tau1(2); end;
if length(tau2) == 1, tau2 = tau2*ones(n-1,1); end

if nargin < 6 || isempty(kappa2), kappa2 = kappa1(2); end;
if length(kappa2) == 1, kappa2 = kappa2*ones(n-1,1); end

% tau1,kappa1,tau2,kappa2

% P = eye(n) - diag(ones(n-1,1),-1);
P = spdiags(ones(n,1)*[-1,1],-1:0,n,n);

% P*diag([tau2(:); 0])*P', tau1
% D = P*diag([tau2(:); 0])*P' + diag(tau1);
D = P*diag(sparse([tau2(:); 0]))*P' + diag(sparse(tau1));

% P*diag([kappa2(:); 0])*P', kappa1
% K = P*diag([kappa2(:); 0])*P' + diag(kappa1);
K = P*diag(sparse([kappa2(:); 0]))*P' + diag(sparse(kappa1));

% M = diag(mu);
M = diag(sparse(mu));

% OUTPUT
coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
