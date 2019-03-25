function [coeffs,fun,F] = orr_sommerfeld(N,R,w)
%ORR_SOMMERFELD   Quartic PEP arising from Orr-Sommerfeld equation.
%  [COEFFS,FUN,F] = nlevp('orr_sommerfeld',N,R,W) returns the coefficients
%  of an N-by-N quartic matrix polynomial
%  P(lambda) = lambda^4*A4 + lambda^3*A3 + lambda^2*A2 + lambda*A1 + A0
%  arising from the spatial stability analysis of the Orr-Sommerfeld equation.
%  The parameter R is the Reynolds number of the problem and w the frequency.
%  The default values are N = 64, R = 5772, W = 0.26943.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2, A3, A4}.
%  FUN is a function handle to evaluate the monomials 1,lambda,..,lambda^4
%  and their derivatives.
%  F is the function handle P(lambda) = lambda^4*A4 + lambda^3*A3 +
%  lambda^2*A2 + lambda*A1 + A0.
%  This problem has the properties pep, parameter-dependent, scalable.

%  Reference:
%  F. Tisseur and N. J. Higham, Structured Pseudospectra for
%  Polynomial Eigenvalue Problems, with Applications,
%  SIAM J. Matrix Anal. Appl., 23(1):187-208, 2001.

if nargin < 3 || isempty(w), w = 0.26943; end
if nargin < 2 || isempty(R), R = 5772; end
if nargin < 1 || isempty(N), N = 64; end
N = N+1;

% 2nd- and 4th-order differentiation matrices:
[D,x] = cheb(N); D2 = D^2; D2 = D2(2:N,2:N);
S = diag([0; 1 ./(1-x(2:N).^2); 0]);
D4 = (diag(1-x.^2)*D^4 - 8*diag(x)*D^3 - 12*D^2)*S;
D4 = D4(2:N,2:N);

% Orr-Sommerfeld operators
i = sqrt(-1);
I = eye(N-1);
coeffs{1} = I;
coeffs{2} = i*R*diag(1-x(2:N).^2);
coeffs{3} = -2*D2-i*w*R*I;
coeffs{4} = -i*R*(diag(1-x(2:N).^2)*D2) - 2i*R*I;
coeffs{5} = D4+i*w*R*D2;

fun = @(lam) nlevp_monomials(lam,4);
F = @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*(coeffs{3}+lam*(coeffs{4} +lam*coeffs{5})));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,x] = cheb(N)
%CHEB     Compute D = differentiation matrix, x = Chebyshev grid.

% From:
% Lloyd N. Trefethen. Spectral Methods in MATLAB. Society for Industrial
% and Applied Mathematics, Philadelphia, PA, USA, 2000.

if N == 0, D = 0; x = 1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));       % off-diagonal entries
D  = D - diag(sum(D,2));                 % diagonal entries
end
