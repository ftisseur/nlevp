function [coeffs,fun,F,xcoeffs] = damped_beam(n)
%DAMPED_BEAM  QEP from simply supported beam damped in the middle.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('damped_beam',N) constructs an N-by-N
%  quadratic matrix polynomial lambda^2*M + lambda*D + K from a finite
%  element model of a beam clamped at both ends with a damper in the middle.
%  N/2 is the number of finite elements.
%  N should be even. Otherwise N-1 is used. Default: N=200.
%  Half of the eigenvalues of the problem are pure imaginary
%  and are eigenvalues of the undamped problem (D = 0).
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  XCOEFFS is the cell {1, L, 1; K, R, M}, where D = L*R' is the 1-rank
%  decomposition of D.
%  This problem has the properties pep, qep, real, symmetric, scalable,
%  sparse, banded (3 bands), low-rank.

%  Reference:
%  N. J. Higham, D. S. Mackey, F. Tisseur and S. D. Garvey.
%  Scaling, sensitivity and stability in the numerical solution of
%  quadratic eigenvalue problems. Internat. J. Numer. Methods Eng.,
%  73(3):344-360, 2008.

if nargin < 1 || isempty(n)
    n = 200;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
            'now targets the true dimension of the problem'])
end

nele = floor(n/2);
n = 2*nele;
% Geometric and material properties of the beam.
width = 0.05;
height = 0.005;
glength = 1;
dlength = glength/nele;
E = 7e10;
I = width*height^3/12;
area = width*height;
rho =  0.674/(area*glength);

damp = 5;

K1 = [12 6*dlength; 6*dlength 4*dlength^2];
K2 = [-12 6*dlength; -6*dlength 2*dlength^2];
K3 = [12 -6*dlength; -6*dlength 4*dlength^2];

M1 = [156 22*dlength; 22*dlength 4*dlength^2];
M2 = [54 -13*dlength ;13*dlength -3*dlength^2];
M3 = [156 -22*dlength; -22*dlength 4*dlength^2];

% u = ones(nele,1);
v1 = [ones(nele,1); 0];
v2 = [0; ones(nele,1)];
Du = spdiags(ones(nele+1,1),1,nele+1,nele+1);
Dv1 = spdiags(v1,0,nele+1,nele+1);
Dv2 = spdiags(v2,0,nele+1,nele+1);

% K = kron(diag(v1),K1)+kron(diag(v2),K3)+kron(diag(u,1),K2)+kron(diag(u,-1),K2');
K = kron(Dv1,K1) + kron(Dv2,K3) + kron(Du,K2) + kron(Du',K2');
ind = [(2:n) n+2];  % Takes care of boundary conditions.
K = (E*I/dlength^3)*K(ind,ind);

% M = kron(diag(v1),M1)+kron(diag(v2),M3)+kron(diag(u,1),M2)+kron(diag(u,-1),M2');
M = kron(Dv1,M1) + kron(Dv2,M3) + kron(Du,M2) + kron(Du',M2');
M = (rho*area*dlength/420)*M(ind,ind);

% D = zeros(n); D(nele,nele) = damp;
D=sparse(nele,nele,damp,n,n);

coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

enele = sparse(zeros(n,1)); 
enele(nele) = 1;
eneler = enele;
eneler(nele) = damp; 
xcoeffs1 = {1,  enele, 1};
xcoeffs2 = {coeffs{1}, eneler', coeffs{3}};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};
end
