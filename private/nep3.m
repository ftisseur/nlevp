function [coeffs, fun, F] = nep3(n, weight, myseed)
%NEP3   NEP with weighted norm coefficients.
%   [COEFFS,FUN,F] = nlevp('nep3', n, weight, myseed) generates an n-by-n
%   nonlinear matrix function F(lambda) = lambda*EYE(n) + A_0 +
%   exp(21*lambda)*A_1 * (lam+4).^0.5*A_2. The matrices A_j are uniformly
%   random distributed with unitary Frobenius norm, except for A_1, whose
%   norm is equal to weight. The default values are n = 10, weight = 1e6,
%   myseed = 42. This problem has the properties nep, scalable,
%   parameter-dependent, random.

%  References: "Created by F. Tisseur, G.M. Negri Porzio"

if nargin < 1 || isempty(n)
    n = 10;
end
if nargin < 2 || isempty(weight)
    weight = 1e6;
end

if nargin < 3 || isempty(myseed)
    myseed = 42;
end

% Save the current rng state
state = rng;
% Set the wanted rng
rng(myseed);

coeffs = cell(1,3);
coeffs{1} = eye(n);
for j = 2:4
coeffs{j} = rand(n)-0.5;
coeffs{j} = coeffs{j}/norm(coeffs{j}, 'fro');
end
coeffs{3} = coeffs{2}*weight;

alpha = 1/2;
fun = @(lam) nep3_fun(lam, alpha);
F = @(lam) lam*coeffs{1} + coeffs{2} + exp(2i*lam)*coeffs{3} + (lam+4).^alpha*coeffs{4};

% Return to the previous rng state
rng(state);
end





function varargout = nep3_fun(lam, alpha)
lam = lam(:);
l = length(lam);

deriv1 = exp(2i*lam);
varargout{1} = [lam, ones(l,1), deriv1, (lam+4).^alpha];
if nargout > 1
    deriv1 = 2i*deriv1;
    derivCoeffs = alpha - [0:nargout-1];
    aux = cumprod(derivCoeffs);
    varargout{2} = [ones(l,1), zeros(l,1), deriv1, aux(1)*(lam*4).^(alpha-1)];
    for j = 3:nargout
        deriv1 = 2i*deriv1;
        varargout{j} = [zeros(l,2), deriv1, aux(j-1)*(lam*4).^(alpha-j+1)];
    end
end

end
