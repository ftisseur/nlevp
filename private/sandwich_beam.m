function [coeffs,fun,F] = sandwich_beam(N,alpha,tau)
%SANDWICH_BEAM  NEP from model of a clamped sandwich beam.
%  [COEFFS,FUN,F] = nlevp('sandwich_beam',N,alpha,tau) returns N-by-N 
%  matrices Ke, M and Kv for the nonlinear eigenvalue problem defined by
%  [Ke-lambda^2*M+g(lambda)*Kv]*x = 0, where
%  g(lambda) = [G0 + Ginf*(i*lambda*tau)^alpha] / [1 + (i*lambda*tau)^alpha],
%  with G0 = 350.4 kPa and Ginf = 3.062 MPa.
%  Only three values of N are allowed: N = 168, N = 840 and N = 3360.
%  The default values are N = 168, alpha = 0.675 and tau = 8.230 ns.
%  The matrices are returned in a cell array: COEFFS = {Ke,M,Kv}.
%  FUN is a function handle to evaluate the functions 1,-lam^2, and g(lambda)
%  and their first derivatives.
%  F is the function handle Ke-lambda^2*M+g(lambda)*Kv.
%  This problem has the properties nep, sparse, banded (6 bands).

%   Reference:
%   R. Van Beeumen, K. Meerbergen and W. Michiels.
%   A rational Krylov method based on Hermite interpolation for nonlinear
%   eigenvalue problems. SIAM J. Sci. Comput., 35(1):A327-A350, 2013.
%
%   Acknowledgment:
%   We are grateful to Natalia Navarrete and Wim Desmet for providing the
%   matrices for this problem.

% define available problem sizes
sizes = [168; 840; 3360];

% set default value for N if not specified
if (nargin < 1) || (isempty(N))
	% take first available size
	N = sizes(1);
else
	% find closest available size
	[~, idx] = min(abs(N - sizes));
	N = sizes(idx);
end

% set default value for alpha if not specified
if (nargin < 2) || (isempty(alpha))
	alpha = .675;
end

% set default value for tau if not specified
if (nargin < 3) || (isempty(tau))
	tau = 8.230e-9;
end

% load matrices
load(sprintf('sandwich_beam_%u.mat', N));


%OUTPUT
coeffs = {Ke, M, Kv};
fun = @(lam) sandwich_beam_fun(lam,alpha,tau);

% Function handle F
G0   = 3.504e5;
Ginf = 3.062e9;
g = @(lam) (G0 + Ginf*(1i*lam*tau).^alpha)./(1 + (1i*lam*tau).^alpha);
F = @(lam) Ke -lam^2*M + g(lam)*Kv;

end


function varargout = sandwich_beam_fun(lam,alpha,tau)

lam = lam(:);
n = length(lam);
G0   = 3.504e5;
Ginf = 3.062e9;
g = (G0 + Ginf*(1i*lam*tau).^alpha)./(1 + (1i*lam*tau).^alpha);

varargout{1} = [ones(n,1),-lam.^2,g];
if nargout >= 2
    g = ((Ginf - G0)*alpha*(1i*lam*tau).^alpha)./(lam.*(1 + (1i*lam*tau).^alpha).^2);
    varargout{2} = [zeros(n,1),-2*lam,g];
end

end
