function [coeffs,fun,F,P] = pdde_stability(n,a0,b0,a1,b1,a2,b2,phi1)
%PDDE_STABILITY  QEP from stability analysis of discretized PDDE.
%  [COEFFS,FUN,F] = nlevp('pdde_stability',N,A0,B0,A1,B1,A2,B2,PHI1)
%  generates the coefficient matrices of a quadratic matrix polynomial
%  lambda^2*E + lambda*F + G of dimension M^2 that is related to the
%  stability analysis of a partial delay-differential equation.
%  M is chosen such that M^2 approximates N.
%  A0, B0, A1, B1, A2, B2, PHI1 are real parameters, with -pi < PHI1 < pi.
%  The default values are N = 225, A0  = 2, B0  = .3, A1  =-2,
%  B1  = .2, A2 = -2, B2 = -.3 and PHI1 = -pi/2.
%  The matrices are returned in a cell array: COEFFS = {G, F, E}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*E + lambda*F + G.
%  [COEFFS,FUN,F,P] = NLEVP('pdde_stability',...) also returns a symmetric
%  permutation matrix P such that E = P*conj(G)*P and F = P*conj(F)*P.
%  This problem has the properties pep, qep, scalable, parameter-dependent,
%  sparse, symmetric, banded (15 bands).

% References
% E. Jarlebring, The Spectrum of Delay-Differential Equations:
%    Numerical Methods, Stability and Perturbation, PhD thesis,
%    TU Braunschweig, Institut Computational Mathematics, Germany, 2008.
% H. Fassbender, N. Mackey, D. S. Mackey and C. Schroeder, Structured
%    Polynomial Eigenproblems Related to Time-Delay Systems, 
%    Electron. Trans. Numer. Anal., 31:306-330, 2008.

if nargin < 1 || isempty(n),  n = 225;    
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
            'now targets the true dimension of the problem']);
end
if nargin < 2 || isempty(a0),   a0   = 2;     end
if nargin < 3 || isempty(b0),   b0   = .3;    end
if nargin < 4 || isempty(a1),   a1   = -2;    end
if nargin < 5 || isempty(b1),   b1   = .2;    end
if nargin < 6 || isempty(a2),   a2   = -2;    end
if nargin < 7 || isempty(b2),   b2   = -.3;   end
if nargin < 8 || isempty(phi1), phi1 = -pi/2; end

m = floor(sqrt(n));
if abs(n-m^2) > abs(n-(m+1)^2)
    m = m+1;
end
n = m^2;

h = pi/(m+1);
x = (1:m)*h;

% Set up matrices for delay equations.
e = ones(m,1);
A0 = spdiags([e, -2*e, e], -1:1, m, m)/(h^2);
A0 = A0 + diag(sparse(a0+b0*sin(x)));
A1 =      diag(sparse(a1+b1*x.*(1-exp(x-pi))));
A2 =      diag(sparse(a2+b2*x.*(pi-x)));

% Convert time delay eqn to quadratic eigenvalue problem
% (lam^2*E + lam*F + P*E*P)*x = 0.
E = kron(speye(m),A2);
gamma = exp(sqrt(-1)*phi1); gamma = gamma/abs(gamma);
F = kron(speye(m),A0-gamma*A1) + kron(A0+gamma*A1,speye(m));

p = reshape(reshape(1:n,m,m)',[],1); % Permutation vector.

% Store coefficients in cell array.
coeffs{3} = E;
coeffs{2} = F;
coeffs{1} = conj(E(p,p));

fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

if nargout > 3
    P = speye(n); P = P(p,:);
end


end
