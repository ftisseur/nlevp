function [coeffs,fun,F,P] = metal_strip
%METAL_STRIP  QEP related to stability of electronic model of metal strip.
%  [COEFFS,FUN,F] = nlevp('metal_strip') generates the coefficient matrices
%  of a quadratic matrix polynomial lambda^2*E + lambda*F + G of
%  dimension 9 that is related to the stability analysis of an electronic
%  model of a metal strip.
%  The matrices are returned in a cell array: COEFFS = {G, F, E}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lambda^2*E + lambda*F + G.
%  [COEFFS,FUN,F,P] = nlevp('metal_strip') also returns a symmetric
%  permutation matrix P such that E = P*conj(G)*P and  F = P*conj(F)*P.
%  This problem has the properties pep, qep, real.

% Reference for model of metal strip:
% A. Bellen, N. Guglielmi and A. E. Ruehli, Methods for Linear Systems
%    of Circuit Delay-Differential Equations of Neutral Type,
%    IEEE Trans. Circuits and Systems - I: Fundamental Theory and Applications,
%    Vol 46 (1), 1999, pp. 212-216
% References for transformation to QEP:
% E. Jarlebring, The Spectrum of Delay-Differential Equations:
%    Numerical Methods, Stability and Perturbation, Phd thesis,
%    TU Braunschweig, Institut Computational Mathematics, Germany, 2008
% H. Fassbender, N. Mackey, D. S. Mackey and C. Schroeder, Structured
%    Polynomial Eigenproblems Related to Time-Delay Systems,
%    Electron. Trans. Numer. Anal., 31:306¡V330, 2008.

% Coefficients for delay differential equation for peec-model of a metal strip.
% D1*x'(t-h)+D0*x'(t)=A1*x(t-h)+A0*x(t)

A0 = 100*[-7,1,2; 3,-9,0; 1,2,-6];
A1 = 100*[1,0,-3; -.5,-.5,-1; -.5,-1.5,0];

D1 = -[-1,5,2; 4,0,3; -2,4,1]/72;
D0 = eye(3);

% Transform time delay system into quadratic eigenvalue problem.
% (lam^2*E + lam*F + P*conj(E)*P)x = 0.
E = kron(conj(D0),A1) + kron(conj(A0),D1);
F = kron(conj(D0),A0) + kron(conj(A0),D0) + kron(D1,A1) + kron(A1,D1);

n = 3;
p = reshape(reshape(1:n^2,n,n)',[],1); % Permutation vector.
% P = eye(n^2); P = P(p,:);            % Permutation matrix.

% Store matrices in cell array.
coeffs{3} = E;
coeffs{2} = F;
coeffs{1} = conj(E(p,p)); % coeffs{1} = P*conj(coeffs{3})*P;

fun = @(lam) nlevp_monomials(lam,2);

F = nlevp_handleQEP(coeffs);

if nargout >= 2
    P = speye(n^2); P = P(p,:);
end


end
