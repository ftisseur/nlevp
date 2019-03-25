function [coeffs, fun, F, sol] = distributed_delay1
%DISTRIBUTED_DELAY1 3-by-3 NEP from distributed delay system.
%  [COEFFS, FUN, F, SOL] = nlevp('distributed_delay1') returns the 
%  3-by-3 matrices I, A1, A2 and A3 in the nonlinear eigenvalue problem
%  R(lambda)x = [-lambda*EYE(3) + A1 + A2*exp(-lambda) +A3*f(lambda)]x = 0.
%  The matrices are returned in a cell array: COEFFS = {-EYE(3),A1,A2,A3}.
%  FUN is a function handle to evaluate the functions lambda,1,
%  exp(-tau*lambda), and f(lambda) and their derivatives.
%  F is the function_handle R(lambda).
%  This problem has the properties nep, real, solution.

%  Reference:
%  E. Jarlebring, W. Michiels and K. Meerbergen.
%  The infinite Arnoldi method and an application to time-delay systems
%  with distributed delays. In Time Delay Systems - Methods, Applications 
%  and New Trends, R. Sipahi, T. Vyhlidal, S.-I. Niculescu, and P. Pepe, 
%  Eds. Number 423 in Lecture Notes in Control and Information Sciences. 
%  Springer-Verlag, Berlin, 229--239, 2012.

% if nargin < 1, quadv_tol = eps; end;

if nargout > 2
    load distributed_delay1; 
end
%   if quadv_tol == eps
%     load distributed_delay1;
%   else
%     error('Solutions only available for QUADV_TOL equal to EPS.');
%   end
% end

A1 = [2.5    2.8   -0.5
      1.8    0.3    0.3
     -2.3   -1.4    3.5];

A2 = [1.7    0.7   -0.3
     -2.4   -2.1   -0.2
      2.0    0.7    0.4];

A3 = [1.4   -1.3    0.4
      1.4    0.7    1.0
      0.6    1.6    1.7];

coeffs = {-eye(3), A1, A2, A3};
fun = @(lam) distributed_delay1_fun(lam);

% Returning the function_handle F
distr_kernel = @(s) exp((s+0.5).^2)-exp(1/4);
f3 = @(lam) integral(@(s) exp(lam*s).*distr_kernel(s),-1,0);
F = @(lam) lam*coeffs{1} + coeffs{2} + coeffs{3}*exp(-lam) + ...
  coeffs{4}*f3(lam);
end

%function varargout = distributed_delay1_fun(lam,quadv_tol)
function varargout = distributed_delay1_fun(lam)
lam = lam(:);
n = length(lam);
distr_kernel = @(s) exp((s+0.5).^2)-exp(1/4);
f3 = integral(@(s) exp(lam*s).*distr_kernel(s),-1,0, 'ArrayValued', true);
varargout{1} = [lam,ones(n,1),exp(-lam),f3];

if nargout >= 2
   f3 = integral(@(s) exp(lam*s).*(s*distr_kernel(s)),-1,0,'ArrayValued', true);
   varargout{2} = [ones(n,1),zeros(n,1),-exp(-lam), f3];
end

for i = 2:nargout-1
   f3 = integral(@(s) exp(lam*s).*((s.^i)*distr_kernel(s)),-1,0,'ArrayValued', true);
   varargout{i+1} = [zeros(n,2), ((-1)^i)*exp(-lam), f3];
end

end
