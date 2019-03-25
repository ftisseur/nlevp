function [coeffs,fun,F,xcoeffs] = photonic_crystal(N, wavevec, eps0, eps1, lorentz)
%PHOTONIC_CRYSTAL  REP from dG-FEM of wave propagation in a periodic medium.
%  [COEFFS,FUN,F,XCOEFFS] = photoniccrystal(N,WAVEVEC,EPS0,EPS1,LORENTZ) generates an
%  N-by-N rational matrix function
%  G-lambda^2*EPS0*M0-lambda^2*eps1(lambda)*M1
%  arising from the dG-finite-element discretization of an eigenvalue problem
%  with periodic boundary conditions describing the propagation of a wave
%  with wave vector WAVEVEC in a 2D, infinite, periodic dielectric medium
%  consisting of two materials. The first material has the (constant) relative
%  permittivity EPS0, whereas the relative permittivity of the second material
%  is given by the Lorentz model
%  eps1(lambda) = EPS1 + sum(lambdaP.^2 ./ (lambda0.^2 - lambda^2 - I*gamma*lambda)),
%  whose coefficients are collected in LORENTZ = [lambdaP.^2; lambda0.^2; gamma].
%  The default values are N = 0, WAVEVEC = [0; 0], EPS0 = 1 for vacuum, and
%  EPS1 = 2, LORENTZ = [[2.5; 1.4; 0.001], [5; 1.6; 0.02]].
%  The coefficient matrices are returned as the cell array COEFFS = {G, M0, M1}.
%  FUN is a function handle to evaluate the corresponding scalar functions 1,
%  -lambda^2*EPS0, and -lambda^2*eps1(lambda) as well as their derivatives.
%  F is the function handle G -lambda^2*EPS0*M0 - lambda^2*eps1(lambda)*M1.
%  XCOEFFS is the cell {1, 1, L1; coeffs{1}, coeffs{2}, R1'} that exploits
%  the low rank of M1.
%  This problem has the properties rep, symmetric, sparse,
%  parameter-dependent, low-rank.

%  References:
%  C. Effenberger, D. Kressner, and C. Engström. Linearization techniques 
%  for band structure calculations in absorbing photonic crystals. 
%  Internat. J. Numer. Methods Engrg., 89(2):180–191, 2012.
%
%  C. Engstrom and M. Wang, M. 2010. Complex dispersion relation
%  calculations with the symmetric interior penalty method. 
%  Internat. J. Numer. Methods Engrg., 84(7):849–863, 2010.

%  Function supplied by Daniel Kressner and Cedric Effenberger.

% define available problem sizes
sizes = [288; 720; 1344; 2160];

% set default value for N if not specified
if (nargin < 1) || (isempty(N))
    % take first available size
	N = sizes(1);
else
	% find closest available size
	[~, idx] = min(abs(N - sizes));
	N = sizes(idx);
end

% set default value for wavevec if not specified
if (nargin < 2) || (isempty(wavevec))
	wavevec = [0; 0];
end

% set default value for eps0 if not specified
if (nargin < 3) || (isempty(eps0))
	eps0 = 1;
end

% set default value for eps1 if not specified
if (nargin < 4) || (isempty(eps1))
	eps1 = 2;
end

% set default value for lorentz if not specified
if (nargin < 5) || (isempty(lorentz))
  lorentz = [2.5 5; 1.4 1.6; 0.001 0.02];
end

% load matrices
load(sprintf('photoniccrystal_%u.mat', N));

% assemble coefficients
coeffs = {ATM + 1i*(wavevec(1)*MDx + wavevec(2)*MDy)' - 1i*(wavevec(1)*MDx + wavevec(2)*MDy) + (wavevec'*wavevec)*(M0+M1), ...
			M0, M1};

% return scalar function
fun = @(lambda) photoniccrystal_fun(lambda, eps0, eps1, lorentz);
% build the function handle F
lambdaPsq = lorentz(1, :);
lambda0sq = lorentz(2, :);
gamm = lorentz(3, :);
eps1fun = @(lam) eps1 + sum(lambdaPsq./(lambda0sq - lam^2 -1i*lam*gamm));
F = @(lam) coeffs{1} - lam^2*eps0*coeffs{2} - lam^2*eps1fun(lam)*coeffs{3};
% building xcoeffs
xcoeffs = {1, 1, L1; coeffs{1}, coeffs{2}, R1'};
end

% helper function to evaluate the scalars
function varargout = photoniccrystal_fun(lambda, eps0, eps1, lorentz)

% convert lambda to column vector
lambda = lambda(:);

% determine number of lambdas supplied
n = length(lambda);

% determne number of Lorentz terms
m = size(lorentz, 2);

% extract coefficients of Lorentz model
lambdaPsq = lorentz(1, :);
lambda0sq = lorentz(2, :);
gamm = lorentz(3, :);

% pre-compute frequently used quantities
lambda2 = lambda.^2;
c = 1i*gamm;
b = lambda*c - lambda0sq(ones(n, 1), :);
a = -1 ./ (b + lambda2(:, ones(m, 1)));

% evaluate scalar functions
varargout{1} = [ones(n, 1), lambda2*(-eps0), ...
	sum(lambdaPsq(ones(n, 1), :) .* (1 + b .* a), 2) + lambda2*(-eps1)];

% evaluate first derivative of scalar functions if requested
if (nargin >= 2)
	
	d = 2*lambda(:, ones(m, 1)) + c(ones(n, 1), :);
	q0 = a;
	q1 = d.*q0.*a;
	varargout{2} = [zeros(n, 1), lambda*(-2*eps0), ...
		sum(lambdaPsq(ones(n, 1), :) .* (b.*q1 + c(ones(n, 1), :).*q0), 2) + lambda*(-2*eps1)];
end

% evaluate second derivative of scalar functions if requested
if (nargout >= 3)
	
	fac = 2;
	q2 = (d.*q1 + q0) .* a;
	q0 = q1;
	q1 = q2;
	varargout{3} = [zeros(n, 1), ones(n, 1)*(-2*eps0), ...
		fac*sum(lambdaPsq(ones(n, 1), :) .* (b.*q1 + c(ones(n, 1), :).*q0), 2) + ones(n, 1)*(-2*eps1)];
end

% evaluate higher derivatives if requested
for k = 3:nargout-1	
	fac = fac * k;
	q2 = (d.*q1 + q0) .* a;
	q0 = q1;
	q1 = q2;
	varargout{k+1} = [zeros(n, 2), ...
		fac*sum(lambdaPsq(ones(n, 1), :) .* (b.*q1 + c(ones(n, 1), :).*q0), 2)];
end

end
