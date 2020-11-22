function [coeffs, fun, F] = time_delay3(n,k,myseed)
%TIME_DELAY3   NEP with high-variance-norm coefficients.
%   [COEFFS,FUN,F] = nlevp('time_delay3', n, k, myseed) generates a n-by-n
%   nonlinear matrix function lambda*EYE(n) +
%   \sum_{j=1}^{k}exp(-lambda*tau_j)A_j, where A_j are random matrices
%   generated with myseed and norm(A_j, inf) = 10^{2j}. The default values
%   are n=10, k=2, myseed=0. The matrices are returned in a cell array of
%   length k+1 COEFFS={EYE(n), A_1, ..., A_k}. FUN is a function handle to
%   evaluate the functions lambda, exp(-lambda*tau_1), ...,
%   exp(-lambda*tau_k) and their derivatives. F is the function handle
%   lambda*EYE(n) + \sum_{j=1}^{k}exp(-lambda*tau_j)A_j. This problem has
%   the properties nep, real, scalable, parameter-dependent, random.

%  Reference: "Created by F. Tisseur"


if nargin < 1 || isempty(n)
    n = 10;
end
 if nargin < 2 || isempty(k)
    k = 2; 
 end
 
 if nargin < 3 || isempty(myseed)
    myseed = 0;
end
 
tau = [1:k];

% Save the current rng state
state = rng;
% Set a standard rng
rng(myseed);

coeffs = cell(1,k+1);
coeffs{1} = eye(n);
for j = 1:k
    aux = rand(n,n);
    aux = aux/norm(aux, inf);
coeffs{j+1} = 10^(2*j)*aux;
end

fun = @(lam) time_delay3_fun(lam, tau);
F = @(lam) time_delay3_F(lam, coeffs, tau);

% Return to the previous rng state
rng(state);
end

function varargout = time_delay3_fun(lam, tau)
lam = lam(:);
l = length(lam);

varargout{1} = [lam exp(-lam*tau)];
if nargout > 1
   deriv = -exp(-lam*tau)*diag(tau);
   varargout{2} = [ones(l,1) deriv];
   for j = 3:nargout
      deriv = -deriv*diag(tau);
       varargout{j} = [zeros(l,1) deriv];
   end
end

end


function F = time_delay3_F(lam, coeffs, tau)
n = size(coeffs{1});
F = lam*eye(n);

% The order of the sum is this one because
% norm(coeffs{j},inf)/norm(coeffs{j+1}, inf) = 1e-2
for j = 1:length(tau)
F = F + exp(-lam*tau(j))*coeffs{j+1};
end

end
