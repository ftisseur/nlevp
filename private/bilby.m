function [coeffs,fun,F] = bilby(beta)
%BILBY   5-by-5 QEP from bilby population model.
%  [COEFFS,FUN,F] = nlevp('bilby',BETA) constructs a 5-by-5 quadratic matrix
%  polynomial lambda^2*A + lambda*B + C arising in a quasi-birth-death
%  process model of the population of the greater bilby,
%  an endangered Australian marsupial.
%  A and C are both singular.  BETA is a parameter (default: 0.5).
%  The matrices are returned in a cell array: COEFFS = {C, B, A}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle C + lambda*B + lambda^2*A.
%  This problem has the properties pep, qep, real, parameter-dependent.

%  Reference:
%  N. J. Higham and H.-M. Kim, Numerical analysis of a quadratic
%  matrix equation, IMA J. Numer. Anal., 20(4): 499-519, 2000.

if nargin < 1, beta = 0.5; end

b = [1, 0.4, 0.25, 0.1, 0];
d = [0, 0.5, 0.55, 0.8, 1];
c = ones(1,5) - b - d;

A0 = mat(0.2,b);
A1 = mat(0.2,c);
A2 = mat(0.2,d);

% Convert to general QME, using \cite[(3.3)]{bblp97}.
C = beta*A0';
B = beta*A1' - eye(5);
A = beta*A2';

coeffs = {C,B,A};

fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
%%%%%%%%%%%%%%%%%%%%%
function A = mat(g,x)
%MAT     Matrix for Bilby model.

x = x(:);
A = zeros(5);
A(:,1) = g*x;
A(5,5) = (1-g)*x(5);
for i=1:4
    A(i,i+1) = (1-g)*x(i);
end
end
